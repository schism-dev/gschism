#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cstdint>
#include <stdexcept>
#include <memory>
#include <string>
#include <cmath>

#include <blend2d.h>

namespace py = pybind11;

static BLFontFace g_font_face;
static BLFont g_font;
static bool g_font_ready = false;

constexpr int textSep = 3;
constexpr BLRgba32 nodeColor(85, 80, 255, 255);
constexpr BLRgba32 sideColor(25, 220, 255, 255);
constexpr BLRgba32 elementColor(20, 250, 26, 255);
constexpr BLRgba32 dryNodeLabelColor(255, 120, 120, 255);
constexpr BLRgba32 drySideLabelColor(255, 175, 95, 255);
constexpr BLRgba32 dryElementLabelColor(255, 90, 170, 255);

static void ensure_text_font() {
    if (g_font_ready) return;

    if (g_font_face.createFromFile("C:/Windows/Fonts/consola.ttf") == BL_SUCCESS) {
        if (g_font.createFromFace(g_font_face, 11.0f) == BL_SUCCESS) {
            g_font_ready = true;
        }
    }
}


void drawSquare(BLContext& ctx, double x, double y) {
    const double size = 2;
    ctx.fillRect(x - size/2, y - size/2, size, size);
}

static BLRgba32 decodePackedColor(double packed, uint8_t alpha = 255u) {
    const auto c = static_cast<uint32_t>(std::llround(packed));
    const uint8_t r = static_cast<uint8_t>((c >> 16) & 0xFFu);
    const uint8_t g = static_cast<uint8_t>((c >> 8) & 0xFFu);
    const uint8_t b = static_cast<uint8_t>(c & 0xFFu);
    return BLRgba32(r, g, b, alpha);
}

static void drawMarker(BLContext& ctx, char marker, double x, double y) {
    constexpr double half = 2.0;
    switch (marker) {
        case '+':
            ctx.strokeLine(BLLine(x - half, y, x + half, y));
            ctx.strokeLine(BLLine(x, y - half, x, y + half));
            break;
        case 'x':
            ctx.strokeLine(BLLine(x - half, y - half, x + half, y + half));
            ctx.strokeLine(BLLine(x - half, y + half, x + half, y - half));
            break;
        case 'o':
            ctx.strokeCircle(BLCircle(x, y, half));
            break;
        case '.':
            drawSquare(ctx, x, y);
            break;
        default:
            // Unsupported marker: do nothing.
            break;
    }
}

// Packed markerCode format (int):
//   bits  0.. 7: marker ASCII (e.g. '+', 'x', 'o', '.')
//   bits  8..15: dash length in pixels (0 => continuous line)
//   bits 16..23: width in half-pixels (2 => 1.0px, 3 => 1.5px, ...)
// Backward compatibility:
//   if markerCode fits in 8 bits, treat it as marker only, dash=0, width=2.
static void decodeMarkerCode(int markerCode, char& marker, int& dashLenPx, int& widthHalfPx) {
    const uint32_t code = static_cast<uint32_t>(markerCode);
    marker = static_cast<char>(code & 0xFFu);

    if ((code & 0xFFFFFF00u) == 0u) {
        dashLenPx = 0;
        widthHalfPx = 2;
        return;
    }

    dashLenPx = static_cast<int>((code >> 8) & 0xFFu);
    widthHalfPx = static_cast<int>((code >> 16) & 0xFFu);
    if (widthHalfPx <= 0) widthHalfPx = 2;
}

static int g_w = 0;
static int g_h = 0;



// frame: uint8 numpy array with shape (h, w, 4), RGBA
void draw_test(py::array_t<std::uint8_t, py::array::c_style> frame) {
    auto info = frame.request();
    if (info.ndim != 3 || info.shape[2] != 4) {
        throw std::runtime_error("Expected frame shape (h, w, 4), dtype=uint8");
    }

}

void drawNode(BLContext& ctx, int nodeIdx, double x, double y, double w, double h, bool drawLabel, bool isDry, bool& drawnNode) {
    if (drawnNode) return;
    drawnNode = true;

    if (!(x >= 0.0 && x < static_cast<double>(w) && y >= 0.0 && y < static_cast<double>(h))) return;

    ctx.setFillStyle(nodeColor); 
    drawSquare(ctx, x, y);
    if (drawLabel) {
        ctx.setFillStyle(isDry ? dryNodeLabelColor : nodeColor);
        ctx.fillUtf8Text(BLPoint(x + textSep, y), g_font, std::to_string(nodeIdx).c_str());
    }
}

void drawMesh(
    py::array_t<std::uint8_t, py::array::c_style> frame,
    py::array_t<double, py::array::c_style> x,
    py::array_t<double, py::array::c_style> y,
    py::array_t<int, py::array::c_style> n1,
    py::array_t<int, py::array::c_style> n2,
    py::array_t<int, py::array::c_style> en,
    py::array_t<int, py::array::c_style> i34,
    py::array_t<int, py::array::c_style> node_isDry,
    py::array_t<int, py::array::c_style> side_isDry,
    py::array_t<int, py::array::c_style> element_isDry,
    py::array_t<double, py::array::c_style> U,
    double zoom,
    double pan_x,
    double pan_y,
    int n_pts,
    py::array_t<double, py::array::c_style> plot_cpu
) {

    auto frameInfo = frame.request();
    if (frameInfo.ndim != 3 || frameInfo.shape[2] != 4) {
        throw std::runtime_error("frame must be uint8[h,w,4]");
    }

    auto xInfo  = x.request();
    auto yInfo  = y.request();
    auto n1Info = n1.request();
    auto n2Info = n2.request();
    auto enInfo = en.request();
    auto i34Info = i34.request();
    auto nodeDryInfo = node_isDry.request();
    auto sideDryInfo = side_isDry.request();
    auto elemDryInfo = element_isDry.request();
    auto UInfo = U.request();

    if (xInfo.ndim != 1 || yInfo.ndim != 1 || n1Info.ndim != 1 || n2Info.ndim != 1) {
        throw std::runtime_error("x, y, n1, n2 must be 1D arrays");
    }
    if (xInfo.shape[0] != yInfo.shape[0]) {
        throw std::runtime_error("x and y must have same length");
    }
    if (n1Info.shape[0] != n2Info.shape[0]) {
        throw std::runtime_error("n1 and n2 must have same length");
    }
    if (nodeDryInfo.ndim != 1 || nodeDryInfo.shape[0] != xInfo.shape[0]) {
        throw std::runtime_error("node_isDry must be 1D and match nodes length");
    }
    if (sideDryInfo.ndim != 1 || sideDryInfo.shape[0] != n1Info.shape[0]) {
        throw std::runtime_error("side_isDry must be 1D and match sides length");
    }
    if (elemDryInfo.ndim != 1 || elemDryInfo.shape[0] != enInfo.shape[1]) {
        throw std::runtime_error("element_isDry must be 1D and match elements length");
    }
    if (i34Info.ndim != 1 || i34Info.shape[0] != enInfo.shape[1]) {
        throw std::runtime_error("i34 must be 1D and match elements length");
    }
    if (zoom <= 0.0) {
        throw std::runtime_error("zoom must be > 0");
    }

    const int h = static_cast<int>(frameInfo.shape[0]);
    const int w = static_cast<int>(frameInfo.shape[1]);
    auto* out = static_cast<uint8_t*>(frameInfo.ptr);
    const intptr_t stride = static_cast<intptr_t>(w * 4);

    const auto* x_ptr  = static_cast<const double*>(xInfo.ptr);
    const auto* y_ptr  = static_cast<const double*>(yInfo.ptr);
    const auto* n1_ptr = static_cast<const int*>(n1Info.ptr);
    const auto* n2_ptr = static_cast<const int*>(n2Info.ptr);
    const auto* en_ptr = static_cast<const int*>(enInfo.ptr);
    const auto* i34_ptr = static_cast<const int*>(i34Info.ptr);
    const auto* nodeDry_ptr = static_cast<const int*>(nodeDryInfo.ptr);
    const auto* sideDry_ptr = static_cast<const int*>(sideDryInfo.ptr);
    const auto* elemDry_ptr = static_cast<const int*>(elemDryInfo.ptr);
    const auto* U_ptr = static_cast<const double*>(UInfo.ptr);
    const py::ssize_t nNodes = xInfo.shape[0];
    const py::ssize_t nSides = n1Info.shape[0];
    const py::ssize_t nElements = enInfo.shape[1];

    static bool *drawnNode = new bool[nNodes];
    for (py::ssize_t i = 0; i < nNodes; ++i) drawnNode[i] = false;


    BLImage img;
    BLResult r = img.createFromData(
        w, h, BL_FORMAT_PRGB32, out, stride,
        BL_DATA_ACCESS_RW, nullptr, nullptr
    );
    if (r != BL_SUCCESS) throw std::runtime_error("Blend2D createFromData failed");

    BLContext ctx(img);
    if (!ctx) throw std::runtime_error("Blend2D context creation failed");

    ctx.setCompOp(BL_COMP_OP_SRC_COPY);
    ctx.setFillStyle(BLRgba32(20, 20, 26, 255));
    ctx.fillAll();


    

    ctx.setCompOp(BL_COMP_OP_SRC_COPY); // opaque lines: faster than SRC_OVER for this case
    ctx.setStrokeStyle(BLRgba32(25, 220, 255, 255));
    ctx.setStrokeWidth(1.0);
    ctx.setStrokeCaps(BL_STROKE_CAP_BUTT);
    ctx.setStrokeJoin(BL_STROKE_JOIN_BEVEL);
    ensure_text_font();
    ctx.setFillStyle(BLRgba32(255, 230, 80, 255));

    auto sx = [zoom, pan_x](double xw) { return xw * zoom + pan_x; };
    auto sy = [zoom, pan_y](double yw) { return yw * zoom + pan_y; };

    auto plotInfo = plot_cpu.request();
    if (plotInfo.ndim != 2 || plotInfo.shape[1] < 7) {
        throw std::runtime_error("plot_cpu must be 2D with 7 columns");
    }
    const auto* plot_ptr = static_cast<const double*>(plotInfo.ptr);
    const py::ssize_t plot_rows = plotInfo.shape[0];
    const py::ssize_t stride0 = plotInfo.strides[0] / sizeof(double);
    const py::ssize_t stride1 = plotInfo.strides[1] / sizeof(double);
    const int n = std::min<int>(n_pts, static_cast<int>(plot_rows));
    BLArray<double> dashPattern;
    for (int i = 0; i < n; ++i) {
        const double x0 = plot_ptr[i * stride0 + 0 * stride1];
        const double y0 = plot_ptr[i * stride0 + 1 * stride1];
        const double x1 = plot_ptr[i * stride0 + 2 * stride1];
        const double y1 = plot_ptr[i * stride0 + 3 * stride1];
        const double packedColor = plot_ptr[i * stride0 + 4 * stride1];
        const int markerCode = static_cast<int>(std::llround(plot_ptr[i * stride0 + 5 * stride1]));
        const int label = static_cast<int>(std::llround(plot_ptr[i * stride0 + 6 * stride1]));

        const double sx0 = sx(x0);
        const double sy0 = sy(y0);
        const double sx1 = sx(x1);
        const double sy1 = sy(y1);

        const BLRgba32 markerColor = decodePackedColor(packedColor);
        const BLRgba32 faintColor = decodePackedColor(packedColor, 110u);
        ctx.setFillStyle(markerColor);
        ctx.setStrokeStyle(markerColor);

        char marker = '.';
        int dashLenPx = 0;
        int widthHalfPx = 2;
        decodeMarkerCode(markerCode, marker, dashLenPx, widthHalfPx);
        ctx.setStrokeWidth(0.5 * static_cast<double>(widthHalfPx));
        if (dashLenPx > 0) {
            dashPattern.clear();
            dashPattern.append(static_cast<double>(dashLenPx));
            dashPattern.append(static_cast<double>(dashLenPx));
            ctx.setStrokeDashOffset(0.0);
            ctx.setStrokeDashArray(dashPattern);
        } else {
            dashPattern.clear();
            ctx.setStrokeDashOffset(0.0);
            ctx.setStrokeDashArray(dashPattern);
        }

        drawMarker(ctx, marker, sx0, sy0);
        drawMarker(ctx, marker, sx1, sy1);

        BLPath seg;
        seg.moveTo(sx0, sy0);
        seg.lineTo(sx1, sy1);
        ctx.strokePath(seg);

        const double dx = sx1 - sx0;
        const double dy = sy1 - sy0;
        const double segLen = std::sqrt(dx * dx + dy * dy);
        const double mx = 0.5 * (sx0 + sx1);
        const double my = 0.5 * (sy0 + sy1);

        // Default direction for degenerate (zero-length) segments.
        double ux = 1.0;
        double uy = 0.0;
        if (segLen > 1e-9) {
            const double invLen = 1.0 / segLen;
            ux = dx * invLen;
            uy = dy * invLen;
        }
        const double px = -uy;
        const double py = ux;

        if (segLen > 1e-9) {
            // Direction chevron (looks like '>') centered at the segment midpoint.
            const double along = 4.0;
            const double arm = 3.0;
            const double tipX = mx + ux * along;
            const double tipY = my + uy * along;
            const double baseX = mx - ux * along;
            const double baseY = my - uy * along;
            ctx.strokeLine(BLLine(baseX + px * arm, baseY + py * arm, tipX, tipY));
            ctx.strokeLine(BLLine(baseX - px * arm, baseY - py * arm, tipX, tipY));
        }

        if (g_font_ready) {
            // Push labels farther as i grows (modulo cycle) and connect to midpoint.
            const double offset = 10.0 + 4.0 * static_cast<double>(i % 11);
            const double tangentShift = 3.0 * static_cast<double>((i % 5) - 2);
            const double lx = mx + px * offset + ux * tangentShift;
            const double ly = my + py * offset + uy * tangentShift;
            ctx.setStrokeStyle(faintColor);
            ctx.setStrokeWidth(1.0);
            dashPattern.clear();
            ctx.setStrokeDashOffset(0.0);
            ctx.setStrokeDashArray(dashPattern);
            ctx.strokeLine(BLLine(mx, my, lx, ly));

            const std::string text = (label != -1) ? std::to_string(label) : std::to_string(i);
            ctx.setFillStyle(markerColor);
            ctx.fillUtf8Text(BLPoint(lx + textSep, ly), g_font, text.c_str());
        }
    }

    auto* base = static_cast<const char*>(enInfo.ptr);
    auto elIdx = [&](py::ssize_t r, py::ssize_t c) -> int {
        return *reinterpret_cast<const int*>(base + r*enInfo.strides[0] + c*enInfo.strides[1]);
    };
    for (py::ssize_t i = 0; i < nElements; ++i) {
        const int i0 = elIdx(0, i);
        const int i1 = elIdx(1, i);
        const int i2 = elIdx(2, i);
        const int i3 = elIdx(3, i);
        const int elem_i34 = i34_ptr[i];
        const bool hasNode4 = (elem_i34 == 4);

        // printf("drawMesh: numElements = %ld, i = %ld, i0 = %d, i1 = %d, i2 = %d, i3 = %d\n", nElements, i, i0, i1, i2, i3);
        
        if (i0 < 0 || i1 < 0 || i2 < 0 ||
                i0 >= nNodes || i1 >= nNodes || i2 >= nNodes) {
                continue;
            }
        if (hasNode4 && (i3 < 0 || i3 >= nNodes)) {
            continue;
        }

        const double e0x = sx(x_ptr[i0]);
        const double e0y = sy(y_ptr[i0]);
        const double e1x = sx(x_ptr[i1]);
        const double e1y = sy(y_ptr[i1]);
        const double e2x = sx(x_ptr[i2]);
        const double e2y = sy(y_ptr[i2]);
        const double e3x = hasNode4 ? sx(x_ptr[i3]) : 0.0;
        const double e3y = hasNode4 ? sy(y_ptr[i3]) : 0.0;

        double elemSize = std::max(std::max(std::abs(e0x - e1x), std::abs(e0y - e1y)),
                                   std::max(std::abs(e0x - e2x), std::abs(e0y - e2y)));
        if (hasNode4) {
            elemSize = std::max(elemSize,
                                std::max(std::max(std::abs(e0x - e3x), std::abs(e0y - e3y)),
                                         std::max(std::abs(e2x - e3x), std::abs(e2y - e3y))));
        }
        const bool drawLabel = elemSize > 50.0;

        drawNode(ctx, i0, e0x, e0y, w, h, drawLabel, nodeDry_ptr[i0] != 0, drawnNode[i0]);
        drawNode(ctx, i1, e1x, e1y, w, h, drawLabel, nodeDry_ptr[i1] != 0, drawnNode[i1]);
        drawNode(ctx, i2, e2x, e2y, w, h, drawLabel, nodeDry_ptr[i2] != 0, drawnNode[i2]);
        if (hasNode4) {
            drawNode(ctx, i3, e3x, e3y, w, h, drawLabel, nodeDry_ptr[i3] != 0, drawnNode[i3]);
        }


        // Draw element centroid
        ctx.setFillStyle(elementColor);

        const double x0 = x_ptr[i0], y0 = y_ptr[i0];
        const double x1 = x_ptr[i1], y1 = y_ptr[i1];
        const double x2 = x_ptr[i2], y2 = y_ptr[i2];
        const double x3 = hasNode4 ? x_ptr[i3] : 0.0;
        const double y3 = hasNode4 ? y_ptr[i3] : 0.0;
        // Element centroid in world coordinates
        const double cxw = hasNode4 ? (x0 + x1 + x2 + x3) / 4.0 : (x0 + x1 + x2) / 3.0;
        const double cyw = hasNode4 ? (y0 + y1 + y2 + y3) / 4.0 : (y0 + y1 + y2) / 3.0;
        // to screen
        const double cx = sx(cxw);
        const double cy = sy(cyw);

        // visible?
        if (cx < 0.0 || cx >= static_cast<double>(w) ||
            cy < 0.0 || cy >= static_cast<double>(h)) {
            continue;
        }
        // draw a small point
        drawSquare(ctx, cx, cy);    
        
        if (drawLabel) {
            const std::string label = std::to_string(static_cast<int>(i));
            ctx.setFillStyle(elemDry_ptr[i] != 0 ? dryElementLabelColor : elementColor);
            ctx.fillUtf8Text(BLPoint(cx + textSep, cy), g_font, label.c_str());
        }
    }

    // Batch all visible side segments into one path, then stroke once.
    static BLPath edges;
    static py::ssize_t edgesReservedSides = -1;
    // Always clear the path each frame, otherwise it will accumulate old lines.
    edges.clear();
    // Only grow capacity when needed (avoid reallocs).
    if (edgesReservedSides != nSides) {
        edges.reserve(static_cast<size_t>(nSides) * 2);
        edgesReservedSides = nSides;
    }

    // Draw sides
    ctx.setFillStyle(sideColor);
    for (py::ssize_t i = 0; i < nSides; ++i) {
        const int i0 = n1_ptr[i];
        const int i1 = n2_ptr[i];


        if (i0 < 0 || i1 < 0 || i0 >= nNodes || i1 >= nNodes) continue;

        const double x0 = x_ptr[i0];
        const double y0 = y_ptr[i0];
        const double x1 = x_ptr[i1];
        const double y1 = y_ptr[i1];

        const double sx0 = sx(x0);
        const double sy0 = sy(y0);
        const double sx1 = sx(x1);
        const double sy1 = sy(y1);

        const double cx = ((sx0 + sx1) / 2.0);
        const double cy = ((sy0 + sy1) / 2.0);

        if ((cx >= 0.0 && cx < static_cast<double>(w) && cy >= 0.0 && cy < static_cast<double>(h))) {
            drawSquare(ctx, cx, cy);
        }

        const double sideSize = std::max(abs(sx0 - sx1), abs(sy0 - sy1));
        const bool drawLabel = sideSize>50.0;
        const bool drawU = sideSize>25.0;

        if (drawLabel &&  (cx >= 0.0 && cx < static_cast<double>(w) &&
                           cy >= 0.0 && cy < static_cast<double>(h))) {
            if (g_font_ready) {
                const std::string label = std::to_string(static_cast<int>(i));
                ctx.setFillStyle(sideDry_ptr[i] != 0 ? drySideLabelColor : sideColor);
                ctx.fillUtf8Text(BLPoint(cx + textSep, cy), g_font, label.c_str());
            }
        }

        if (drawU) {
            const double ux = U_ptr[i*2];
            const double uy = U_ptr[i*2+1];
            ctx.setStrokeStyle(BLRgba32(0, 255, 0, 181));
            ctx.setStrokeWidth(1.0);
            if (std::abs(ux) < 10.0 && std::abs(uy) < 10.0) {
                ctx.strokeLine(BLLine(cx, cy, cx + 10 * ux, cy + 10 * uy));
            }
        }

        // ctx.strokeLine(BLLine(sx(x0), sy(y0), sx(x1), sy(y1)));

        if ((sx0 >= 0.0 && sx0 < static_cast<double>(w) && sy0 >= 0.0 && sy0 < static_cast<double>(h)) ||
            (sx1 >= 0.0 && sx1 < static_cast<double>(w) && sy1 >= 0.0 && sy1 < static_cast<double>(h)) ||
            (cx >= 0.0 && cx < static_cast<double>(w) && cy >= 0.0 && cy < static_cast<double>(h))) {

            edges.moveTo(sx0, sy0);
            edges.lineTo(sx1, sy1);
        }
    }
    ctx.setStrokeStyle(sideColor);                  // sides
    ctx.setStrokeWidth(1.0);
    ctx.strokePath(edges);


    ctx.end();
 
}

PYBIND11_MODULE(render_cpp, m) {
    m.doc() = "canvas_ity renderer";
    m.def("draw_test", &draw_test, py::arg("frame"));
    m.def(
        "drawMesh",
        &drawMesh,
        py::arg("frame"),
        py::arg("x"),
        py::arg("y"),
        py::arg("n1"),
        py::arg("n2"),
        py::arg("en"),
        py::arg("i34"),
        py::arg("node_isDry"),
        py::arg("side_isDry"),
        py::arg("element_isDry"),
        py::arg("U"),
        py::arg("zoom") ,
        py::arg("pan_x"),
        py::arg("pan_y"),
        py::arg("n_pts"),
        py::arg("plot_cpu")
    );
}