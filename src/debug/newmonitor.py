from __future__ import annotations
import glob
import importlib.util
import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11
from python.schism_consts import VISTA
import sys
import time
import numpy as np
import cupy as cp
from pathlib import Path
import debug.monitor as mon

def _node_pointer_info():
    """Return (node_base_ptr, node_stride_bytes) for pointer->index conversion."""
    var_code = list(mon.varInfoDict.keys()).index("x")  # any node var
    offset = cp.zeros(1, dtype=cp.uint64)
    addr = cp.zeros(1, dtype=cp.uint64)
    size_node = cp.zeros(1, dtype=cp.uint64)
    size_type = cp.zeros(1, dtype=cp.uint64)
    mon.cufunc.getNodeVarAddr(var_code, addr, offset, size_node, size_type)
    return int(addr.item()), int(size_node.item())


def _get_var_with_fallback(getter, names: tuple[str, ...]):
    """Try multiple variable names and return the first available view."""
    last_exc = None
    for name in names:
        try:
            return getter(name)
        except Exception as exc:
            last_exc = exc
    raise RuntimeError(f"None of these variables are available: {names}") from last_exc


def extract_nodes():
    node_xy = mon.getNodeVar("x").get()   # (2, numNodes)
    x = np.ascontiguousarray(node_xy[0, :], dtype=np.float64)
    y = np.ascontiguousarray(node_xy[1, :], dtype=np.float64)
    node_dry = _get_var_with_fallback(mon.getNodeVar, ("isdry", "isDry")).get()
    node_dry = np.ascontiguousarray(node_dry.astype(np.int32), dtype=np.int32)
    return {"x": x, "y": y, "isDry": node_dry}


def extract_sides(n_nodes: int, node_base: int, node_stride: int):
    n1_ptr = mon.getSideVar("n1").astype(cp.uint64).get()
    n2_ptr = mon.getSideVar("n2").astype(cp.uint64).get()
    n1 = ((n1_ptr - node_base) // node_stride).astype(np.int64)
    n2 = ((n2_ptr - node_base) // node_stride).astype(np.int64)

    side_dry_all = _get_var_with_fallback(mon.getSideVar, ("isdry", "isDry")).get().astype(np.int32)
    side_dry_all = np.asarray(side_dry_all, dtype=np.int32)

    valid = (n1 >= 0) & (n1 < n_nodes) & (n2 >= 0) & (n2 < n_nodes)
    n1 = np.ascontiguousarray(n1[valid], dtype=np.int32)
    n2 = np.ascontiguousarray(n2[valid], dtype=np.int32)
    side_dry = np.ascontiguousarray(side_dry_all[valid], dtype=np.int32)
    return {"n1": n1, "n2": n2, "isDry": side_dry}


def extract_elements(n_nodes: int, node_base: int, node_stride: int):
    # shape (4, numElements), uint64 pointers to Node*
    face_ptrs = mon.getElementVar("FaceIndices").astype(cp.uint64).get()
    elem_nodes = ((face_ptrs - node_base) // node_stride).astype(np.int32)
    i0 = elem_nodes[0, :]
    i1 = elem_nodes[1, :]
    i2 = elem_nodes[2, :]
    i3 = elem_nodes[3, :]

    base_valid = (
        (i0 >= 0) & (i0 < n_nodes) &
        (i1 >= 0) & (i1 < n_nodes) &
        (i2 >= 0) & (i2 < n_nodes)
    )
    has_node4 = (
        base_valid &
        (i3 >= 0) & (i3 < n_nodes) &
        (i3 != i0) & (i3 != i1) & (i3 != i2)
    )
    # For triangles, nodes[3] is not guaranteed to be initialized.
    elem_nodes[3, ~has_node4] = -1
    elem_i34 = np.where(has_node4, 4, 3).astype(np.int32)

    elem_dry = _get_var_with_fallback(mon.getElementVar, ("isdry", "isDry")).get().astype(np.int32)
    elem_dry = np.asarray(elem_dry, dtype=np.int32)

    # Normalize element dry array to shape (ne,)
    if elem_dry.ndim > 1:
        elem_dry = elem_dry.reshape(-1)

    return {
        "en": np.ascontiguousarray(elem_nodes, dtype=np.int32),
        "i34": np.ascontiguousarray(elem_i34, dtype=np.int32),
        "isDry": np.ascontiguousarray(elem_dry, dtype=np.int32),
    }


def extract_mesh_data():
    node_base, node_stride = _node_pointer_info()
    nodes = extract_nodes()
    sides = extract_sides(nodes["x"].shape[0], node_base, node_stride)
    elements = extract_elements(nodes["x"].shape[0], node_base, node_stride)
    return {"nodes": nodes, "sides": sides, "elements": elements}


def compile_and_import(module_name: str, cpp_file: str):
    cpp_file = os.path.abspath(cpp_file)
    if not os.path.exists(cpp_file):
        raise FileNotFoundError(f"C++ source file not found: {cpp_file}")
    ext_modules = [
        Extension(
            module_name,
            sources=[cpp_file],
            include_dirs=[ pybind11.get_include(),
                pybind11.get_include(user=True),
                r"C:\Users\jmgonzalez\apps\vcpkg\installed\x64-windows\include",
            ],
            library_dirs=[
                r"C:\Users\jmgonzalez\apps\vcpkg\installed\x64-windows\lib",
            ],
            libraries=["blend2d"],
            language="c++",
            extra_compile_args=["/std:c++17", "/Zi", "/O2", "/EHsc"] if os.name == "nt" else ["-std=c++17", "-O3", "-fPIC"],
        )
    ]
    class BuildExt(build_ext):
        def build_extensions(self):
            super().build_extensions()
    setup(
        name=module_name,
        ext_modules=ext_modules,
        cmdclass={"build_ext": BuildExt},
        script_args=["build_ext", "--inplace"],
    )
    pattern = f"{module_name}*.pyd" if os.name == "nt" else f"{module_name}*.so"
    module_files = glob.glob(pattern)
    if not module_files:
        raise ImportError(f"Compiled shared library for {module_name} not found.")
    module_path = max(module_files, key=os.path.getmtime)

    def _resolve_vcpkg_bin():
        # 1) explicit override
        root = os.environ.get("VCPKG_ROOT")
        if root:
            cand = os.path.join(root, "installed", "x64-windows", "bin")
            if os.path.isdir(cand):
                return cand

        # 2) common local fallback relative to user profile
        cand = os.path.join(os.path.expanduser("~"), "apps", "vcpkg", "installed", "x64-windows", "bin")
        if os.path.isdir(cand):
            return cand

        return None

    vcpkg_bin = _resolve_vcpkg_bin()
    if not vcpkg_bin:
        raise RuntimeError("Could not find vcpkg bin. Set VCPKG_ROOT or install vcpkg in ~/apps/vcpkg.")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(vcpkg_bin)
    os.environ["PATH"] = vcpkg_bin + os.pathsep + os.environ.get("PATH", "")

    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if not VISTA:
    from PySide6 import QtWidgets, QtGui, QtCore

    class MeshWidget(QtWidgets.QWidget):
        def __init__(self, module, parent=None):
            super().__init__(parent)
            self.module = module
            self.frame = None
            self.last_ms = 0.0
            self.zoom = 10.0
            self.pan_x = 0.0
            self.pan_y = 0.0
            self._panning = False
            self._last_mouse_pos = QtCore.QPointF()
            self._mouse_screen_pos = QtCore.QPointF(0.0, 0.0)
            self._mouse_world_pos = (0.0, 0.0)
            self.setMinimumSize(300, 200)
            self.setMouseTracking(True)

            self.mesh_ready = False
            self.x = self.y = self.n1 = self.n2 = None
            self.side_cx = None
            self.side_cy = None
            self.side_ids = None
            self.label_font = QtGui.QFont("Consolas", 7)
            self._view_initialized = False
            self._fit_zoom = 1.0
            self._focus_zoom = 1.0
            self._init_mesh_once()


        def _init_mesh_once(self):
            if self.mesh_ready:
                return
            mesh = extract_mesh_data()
            nodes = mesh["nodes"]
            sides = mesh["sides"]
            elements = mesh["elements"]

            self.x = np.ascontiguousarray(nodes["x"], dtype=np.float64)
            self.y = np.ascontiguousarray(nodes["y"], dtype=np.float64)
            self.node_isDry = np.ascontiguousarray(nodes["isDry"], dtype=np.int32)

            # Keep only valid connectivity and preserve original side ids for labels.
            n1 = np.ascontiguousarray(sides["n1"], dtype=np.int32)
            n2 = np.ascontiguousarray(sides["n2"], dtype=np.int32)
            side_ids = np.arange(n1.shape[0], dtype=np.int32)
            valid = (n1 >= 0) & (n2 >= 0) & (n1 < self.x.shape[0]) & (n2 < self.x.shape[0])
            n1 = n1[valid]
            n2 = n2[valid]
            side_ids = side_ids[valid]

            self.n1 = np.ascontiguousarray(n1, dtype=np.int32)
            self.n2 = np.ascontiguousarray(n2, dtype=np.int32)
            self.side_isDry = np.ascontiguousarray(sides["isDry"][valid], dtype=np.int32)
            self.side_ids = np.ascontiguousarray(side_ids, dtype=np.int32)

            self.en = np.ascontiguousarray(elements["en"], dtype=np.int32)
            self.i34 = np.ascontiguousarray(elements["i34"], dtype=np.int32)
            self.element_isDry = np.ascontiguousarray(elements["isDry"], dtype=np.int32)

            self.mesh_ready = True

        def _fit_view_to_mesh(self):
            if not self.mesh_ready or self.x is None or self.y is None:
                return
            if self.x.size == 0 or self.y.size == 0:
                return

            w = max(1, self.width())
            h = max(1, self.height())

            xmin = float(np.min(self.x))
            xmax = float(np.max(self.x))
            ymin = float(np.min(self.y))
            ymax = float(np.max(self.y))
            mesh_w = max(1e-12, xmax - xmin)
            mesh_h = max(1e-12, ymax - ymin)

            margin = 0.95  # keep a small border around the mesh
            self.zoom = 5*margin * min(w / mesh_w, h / mesh_h)
            self.zoom = 1.0

            cx = 0.5 * (xmin + xmax)
            cy = 0.5 * (ymin + ymax)
            self.pan_x = 0.5 * w - cx * self.zoom
            self.pan_y = 0.5 * h - cy * self.zoom
            self._fit_zoom = self.zoom
            self._focus_zoom = self.zoom * 8.0
            self._view_initialized = True

            self.U = mon.getSideVar("u")[:,-1,:].get()
            self.U = np.ascontiguousarray(self.U, dtype=np.float64)

        def _world_to_screen(self, wx: float, wy: float):
            return wx * self.zoom + self.pan_x, wy * self.zoom + self.pan_y

        def _screen_to_world(self, sx: float, sy: float):
            return (sx - self.pan_x) / self.zoom, (sy - self.pan_y) / self.zoom

        def _update_mouse_world(self, pos: QtCore.QPointF):
            self._mouse_screen_pos = pos
            self._mouse_world_pos = self._screen_to_world(pos.x(), pos.y())

        def zoom_to(self, target_kind: str, target_index: int) -> bool:
            if not self.mesh_ready:
                return False
            if target_index < 0:
                return False

            wx = wy = None
            n_nodes = self.x.shape[0]
            n_sides = self.n1.shape[0]
            n_elems = self.en.shape[1] if self.en is not None and self.en.ndim == 2 else 0

            if target_kind == "node":
                if target_index >= n_nodes:
                    return False
                wx = float(self.x[target_index])
                wy = float(self.y[target_index])
            elif target_kind == "side":
                if target_index >= n_sides:
                    return False
                i0 = int(self.n1[target_index])
                i1 = int(self.n2[target_index])
                wx = 0.5 * (float(self.x[i0]) + float(self.x[i1]))
                wy = 0.5 * (float(self.y[i0]) + float(self.y[i1]))
            elif target_kind == "element":
                if target_index >= n_elems:
                    return False
                nodes = self.en[:, target_index].astype(np.int64)
                valid = (nodes >= 0) & (nodes < n_nodes)
                if not np.any(valid):
                    return False
                nodes = nodes[valid]
                wx = float(np.mean(self.x[nodes]))
                wy = float(np.mean(self.y[nodes]))
            else:
                return False

            w = max(1, self.width())
            h = max(1, self.height())
            self.zoom = max(self.zoom, self._focus_zoom, self._fit_zoom)
            self.pan_x = 0.5 * w - wx * self.zoom
            self.pan_y = 0.5 * h - wy * self.zoom
            self.update()
            return True

        def _ensure_frame(self):
            w = max(1, self.width())
            h = max(1, self.height())
            if self.frame is None or self.frame.shape[1] != w or self.frame.shape[0] != h:
                self.frame = np.zeros((h, w, 4), dtype=np.uint8)

        def resizeEvent(self, event):
            self._ensure_frame()
            if not self._view_initialized:
                self._fit_view_to_mesh()
            self.update()  # request redraw at new size
            super().resizeEvent(event)

        def paintEvent(self, event):
            self._ensure_frame()

            n_pts = int(mon.debugPlotIdx.get()[0])          # or int(mon.debugPlotIdx.item())
            plot_cpu = cp.asnumpy(mon.debugPlot[:n_pts, :]) # shape (n_pts, 7), np.float64



            t0 = time.perf_counter()
            self.module.drawMesh(
                self.frame,
                self.x,
                self.y,
                self.n1,
                self.n2,
                self.en,
                self.i34,
                self.node_isDry,
                self.side_isDry,
                self.element_isDry,
                self.U,
                self.zoom,
                self.pan_x,
                self.pan_y,
                n_pts,
                plot_cpu,
            )
            t1 = time.perf_counter()
            self.last_ms = (t1 - t0) * 1000.0

            h, w, _ = self.frame.shape
            img = QtGui.QImage(self.frame.data, w, h, 4 * w, QtGui.QImage.Format.Format_RGBA8888)

            painter = QtGui.QPainter(self)
            painter.drawImage(0, 0, img)
            painter.setPen(QtGui.QColor(255, 255, 255))
            painter.drawText(
                10, 20,
                f"{self.last_ms:.1f} ms | zoom={self.zoom:.3f} pan=({self.pan_x:.1f},{self.pan_y:.1f})"
            )
            wx, wy = self._mouse_world_pos
            painter.drawText(
                10, 38,
                f"mouse(mesh)=({wx:.3f},{wy:.3f}) screen=({self._mouse_screen_pos.x():.1f},{self._mouse_screen_pos.y():.1f})"
            )

            # Draw side index labels only when zoomed in enough and visible on screen.
            if self.zoom > 0.05 and self.side_cx is not None and self.side_cy is not None:
                scale = self.zoom / 1000.0
                sx = self.side_cx * scale + self.pan_x
                sy = self.side_cy * scale + self.pan_y
                mask = (sx >= 0.0) & (sx < w) & (sy >= 0.0) & (sy < h)
                if np.any(mask):
                    painter.setFont(self.label_font)
                    painter.setPen(QtGui.QColor(255, 230, 80))
                    sx_vis = sx[mask]
                    sy_vis = sy[mask]
                    ids_vis = self.side_ids[mask]
                    for px, py, sid in zip(sx_vis, sy_vis, ids_vis):
                        painter.drawText(int(px), int(py), str(int(sid)))
            painter.end()

        def mousePressEvent(self, event):
            if event.button() == QtCore.Qt.MouseButton.LeftButton:
                self._panning = True
                self._last_mouse_pos = event.position()
                self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)
                event.accept()
                return
            super().mousePressEvent(event)

        def mouseMoveEvent(self, event):
            self._update_mouse_world(event.position())
            if self._panning:
                pos = event.position()
                delta = pos - self._last_mouse_pos
                self._last_mouse_pos = pos
                self.pan_x += delta.x()
                self.pan_y += delta.y()
                self.update()
                event.accept()
                return
            super().mouseMoveEvent(event)

        def mouseReleaseEvent(self, event):
            if event.button() == QtCore.Qt.MouseButton.LeftButton and self._panning:
                self._panning = False
                self.setCursor(QtCore.Qt.CursorShape.ArrowCursor)
                event.accept()
                return
            super().mouseReleaseEvent(event)

        def wheelEvent(self, event):
            self._update_mouse_world(event.position())
            angle = event.angleDelta().y()
            if angle == 0:
                event.ignore()
                return

            # Zoom around the cursor position (screen-space pivot).
            zoom_factor = 1.25 if angle > 0 else (1.0 / 1.25)
            old_zoom = self.zoom
            new_zoom = old_zoom * zoom_factor
            if new_zoom <= 0.0 or not np.isfinite(new_zoom):
                event.accept()
                return

            pivot = event.position()
            wx = (pivot.x() - self.pan_x) / old_zoom
            wy = (pivot.y() - self.pan_y) / old_zoom
            self.zoom = new_zoom
            self.pan_x = pivot.x() - wx * new_zoom
            self.pan_y = pivot.y() - wy * new_zoom
            self._update_mouse_world(pivot)
            self.update()
            event.accept()


    def show_window(module):
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

        dialog = QtWidgets.QDialog()
        dialog.setWindowTitle("render_cpp output")
        dialog.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
        layout = QtWidgets.QVBoxLayout(dialog)

        canvas = MeshWidget(module)
        layout.addWidget(canvas, 1)

        controls_widget = QtWidgets.QWidget()
        controls_widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Fixed)
        controls = QtWidgets.QHBoxLayout(controls_widget)
        controls.setContentsMargins(0, 0, 0, 0)
        target_combo = QtWidgets.QComboBox()
        target_combo.addItems(["side", "node", "element"])
        target_index = QtWidgets.QLineEdit()
        target_index.setPlaceholderText("index")
        target_index.setText("1276362")
        target_index.setValidator(QtGui.QIntValidator(0, 2_000_000_000, target_index))
        zoom_btn = QtWidgets.QPushButton("Zoom To")
        status = QtWidgets.QLabel("")

        def on_zoom_to():
            text = target_index.text().strip()
            if not text:
                status.setText("Enter an index")
                return
            ok = canvas.zoom_to(target_combo.currentText(), int(text))
            status.setText("" if ok else "Invalid target/index")

        zoom_btn.clicked.connect(on_zoom_to)
        target_index.returnPressed.connect(on_zoom_to)

        controls.addWidget(target_combo)
        controls.addWidget(target_index)
        controls.addWidget(zoom_btn)
        controls.addWidget(status, 1)
        layout.addWidget(controls_widget, 0)
        layout.setStretch(0, 1)
        layout.setStretch(1, 0)

        dialog.resize(900, 600)
        dialog.show()
        app.exec()


module = None


def show():
    """Open the mesh debug window (desktop only; no-op on VISTA cluster)."""
    if VISTA:
        return

    global module

    if module is None:
        module = compile_and_import("render_cpp", str(Path(__file__).with_name("render.cpp")))
    show_window(module)