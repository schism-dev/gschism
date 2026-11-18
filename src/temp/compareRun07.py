"""
Compare SCHISM (CPU) and gschism (GPU) outputs of a test-suite run at 10 points.

    python compareRun07.py [RUN]      RUN = RUN07a (default), RUN06a, RUN06b, RUN06c, ...

Reads out2d (elevation), horizontalVelX/Y, verticalVelocity and salinity,
snaps 10 points spread over the domain to the nearest grid nodes, and writes:

    compare<RUN>_map.png        point locations over the bathymetry
    compare<RUN>_eta.png        elevation time series
    compare<RUN>_velocity.png   u, v profiles at PROFILE_HOURS (rows) x points (columns)
    compare<RUN>_wvel.png       w profiles, same layout
    compare<RUN>_salinity.png   salinity profiles, same layout (skipped if either run
                                has no salinity_1.nc)

CPU  : .../<RUN>/outputs
GPU  : .../<RUN>/gpuoutput
"""

import os
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

RUN = sys.argv[1] if len(sys.argv) > 1 else "RUN07a"
SUITE_DIR = r"C:\Users\jmgonzalez\src\schism_test_suit"
CPU_DIR = os.path.join(SUITE_DIR, RUN, "outputs")
GPU_DIR = os.path.join(SUITE_DIR, RUN, "gpuoutput")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
PREFIX = "compare" + RUN

# 5 x 2 points, inset from the 0..64000 x 0..20000 domain.
POINT_XY = [(x, y) for y in (5000.0, 15000.0) for x in (8000.0, 20000.0, 32000.0, 44000.0, 56000.0)]

FILL_ABS = 1.0e20

# Profile figures use the records closest to these times.
PROFILE_HOURS = (3.0, 6.0, 9.0, 12.0)


def openVar(folder, filename, varname):
    ds = nc.Dataset(os.path.join(folder, filename))
    return ds, ds.variables[varname]


def nearestNodes(x, y, points):
    idx = []
    for px, py in points:
        d2 = (x - px)**2 + (y - py)**2
        idx.append(int(np.argmin(d2)))
    return idx


def wetMask(z, bottomIndex):
    # bottom_index_node is 1-based. Layers below it are under the bed.
    k = np.arange(z.shape[-1])
    wet = k[None, :] >= (bottomIndex[:, None] - 1)
    wet &= np.abs(z) < FILL_ABS
    return wet


def main():
    cpu2d, elevC = openVar(CPU_DIR, "out2d_1.nc", "elevation")
    gpu2d, elevG = openVar(GPU_DIR, "out2d_1.nc", "elevation")

    x = np.array(cpu2d["SCHISM_hgrid_node_x"][:])
    y = np.array(cpu2d["SCHISM_hgrid_node_y"][:])
    xg = np.array(gpu2d["SCHISM_hgrid_node_x"][:])
    yg = np.array(gpu2d["SCHISM_hgrid_node_y"][:])
    if x.shape != xg.shape or np.max(np.abs(x - xg)) > 1.0e-3 or np.max(np.abs(y - yg)) > 1.0e-3:
        raise RuntimeError("CPU and GPU node coordinates do not match")

    timeC = np.array(cpu2d["time"][:])
    timeG = np.array(gpu2d["time"][:])
    nTime = min(len(timeC), len(timeG))
    hours = timeC[:nTime] / 3600.0

    nodes = nearestNodes(x, y, POINT_XY)
    print("node comparisons (target -> node, x, y):")
    for (px, py), n in zip(POINT_XY, nodes):
        print(f"  ({px:8.0f}, {py:8.0f}) -> node {n:5d}  ({x[n]:10.2f}, {y[n]:10.2f})")

    # Point map over the bathymetry (Delaunay of the nodes is fine for this rectangular domain).
    depth = np.array(cpu2d["depth"][:])
    fig, a = plt.subplots(figsize=(12, 4.5))
    cs = a.tricontourf(x / 1000.0, y / 1000.0, depth, levels=20, cmap="Blues")
    fig.colorbar(cs, ax=a, label="depth (m)")
    a.plot(x[nodes] / 1000.0, y[nodes] / 1000.0, "r^", ms=8)
    for n in nodes:
        a.annotate(f"{n}", (x[n] / 1000.0, y[n] / 1000.0), xytext=(5, 5),
                   textcoords="offset points", fontsize=9, color="r")
    a.set_aspect("equal")
    a.set_xlabel("x (km)")
    a.set_ylabel("y (km)")
    a.set_title(RUN + " comparison points (node numbers)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, PREFIX + "_map.png"), dpi=120)

    etaC = np.array(elevC[:nTime, nodes])
    etaG = np.array(elevG[:nTime, nodes])

    fig, axes = plt.subplots(2, 5, figsize=(16, 6), sharex=True)
    for a, n, (px, py), ec, eg in zip(axes.ravel(), nodes, POINT_XY, etaC.T, etaG.T):
        a.plot(hours, ec, "k-", lw=1.2, label="CPU")
        a.plot(hours, eg, "r--", lw=1.0, label="GPU")
        a.set_title(f"node {n}\n({x[n]/1000:.1f}, {y[n]/1000:.1f}) km", fontsize=9)
        a.grid(True, alpha=0.3)
        a.text(0.02, 0.95, f"max|d|={np.max(np.abs(ec-eg)):.3e}", transform=a.transAxes,
               va="top", fontsize=8)
    axes[0, 0].legend(fontsize=8)
    for a in axes[-1, :]:
        a.set_xlabel("hours")
    for a in axes[:, 0]:
        a.set_ylabel("elevation (m)")
    fig.suptitle(RUN + " elevation: CPU (outputs) vs GPU (gpuoutput)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, PREFIX + "_eta.png"), dpi=120)

    # Profiles at the records closest to PROFILE_HOURS. Each run keeps its own z.
    times = [int(np.argmin(np.abs(hours - h))) for h in PROFILE_HOURS]
    cpuZ, zCvar = openVar(CPU_DIR, "zCoordinates_1.nc", "zCoordinates")
    gpuZ, zGvar = openVar(GPU_DIR, "zCoordinates_1.nc", "zCoordinates")
    cpuU, uCvar = openVar(CPU_DIR, "horizontalVelX_1.nc", "horizontalVelX")
    gpuU, uGvar = openVar(GPU_DIR, "horizontalVelX_1.nc", "horizontalVelX")
    cpuV, vCvar = openVar(CPU_DIR, "horizontalVelY_1.nc", "horizontalVelY")
    gpuV, vGvar = openVar(GPU_DIR, "horizontalVelY_1.nc", "horizontalVelY")
    cpuW, wCvar = openVar(CPU_DIR, "verticalVelocity_1.nc", "verticalVelocity")
    gpuW, wGvar = openVar(GPU_DIR, "verticalVelocity_1.nc", "verticalVelocity")
    hasSalt = all(os.path.exists(os.path.join(d, "salinity_1.nc")) for d in (CPU_DIR, GPU_DIR))
    if hasSalt:
        cpuS, sCvar = openVar(CPU_DIR, "salinity_1.nc", "salinity")
        gpuS, sGvar = openVar(GPU_DIR, "salinity_1.nc", "salinity")
    else:
        print("salinity_1.nc missing in CPU or GPU output: salinity skipped")

    botC = np.array(cpu2d["bottom_index_node"][:])[nodes]
    botG = np.array(gpu2d["bottom_index_node"][:])[nodes]

    # Arrays are (nTimes, nPoints, nLayers).
    def readAt(var):
        return np.array([np.array(var[it, nodes, :]) for it in times])

    zC, zG = readAt(zCvar), readAt(zGvar)
    uC, uG = readAt(uCvar), readAt(uGvar)
    vC, vG = readAt(vCvar), readAt(vGvar)
    wC, wG = readAt(wCvar), readAt(wGvar)
    if hasSalt:
        sC, sG = readAt(sCvar), readAt(sGvar)
    wetC = np.array([wetMask(z, botC) for z in zC])
    wetG = np.array([wetMask(z, botG) for z in zG])

    print(f"\neta  max |CPU-GPU| over points and times: {np.max(np.abs(etaC-etaG)):.4e} m")
    print("profiles at t =", ", ".join(f"{hours[it]:.2f} h (record {it})" for it in times))

    # series: (cpuData, gpuData, cpuStyle, gpuStyle, cpuLabel, gpuLabel)
    # All panels share one x range: min/max of the wet CPU and GPU values over every
    # point and time, padded 5%. minXMax: lower bound for the right x limit, so
    # near-zero profiles are not blown up.
    def profileFig(series, xlabel, title, filename, minXMax=None):
        xMin, xMax = np.inf, -np.inf
        for dC, dG, *rest in series:
            xMin = min(xMin, dC[wetC].min(), dG[wetG].min())
            xMax = max(xMax, dC[wetC].max(), dG[wetG].max())
        if minXMax is not None:
            xMax = max(xMax, minXMax)
        pad = 0.05*(xMax - xMin) if xMax > xMin else 0.05*max(abs(xMax), 1.0e-12)
        xLim = (xMin - pad, xMax + pad)

        fig, axes = plt.subplots(len(times), len(nodes), figsize=(26, 3.0*len(times)),
                                 sharex=True, sharey=True, squeeze=False)
        for r, it in enumerate(times):
            for i, n in enumerate(nodes):
                a = axes[r, i]
                wc, wg = wetC[r, i], wetG[r, i]
                maxDiff = 0.0
                for dC, dG, stC, stG, lbC, lbG in series:
                    a.plot(dC[r, i][wc], zC[r, i][wc], label=lbC, **stC)
                    a.plot(dG[r, i][wg], zG[r, i][wg], label=lbG, **stG)
                    if np.array_equal(wc, wg):
                        maxDiff = max(maxDiff, np.max(np.abs(dC[r, i][wc] - dG[r, i][wg])))
                a.text(0.02, 0.98, f"t = {hours[it]:.2f} h\nmax|d|={maxDiff:.2e}",
                       transform=a.transAxes, va="top", fontsize=7)
                a.set_xlim(xLim)
                a.grid(True, alpha=0.3)
                a.tick_params(labelsize=7)
                if r == 0:
                    a.set_title(f"node {n}\n({x[n]/1000:.0f}, {y[n]/1000:.0f}) km", fontsize=8)
                if r == len(times) - 1:
                    a.set_xlabel(xlabel, fontsize=8)
            axes[r, 0].set_ylabel(f"t = {hours[it]:.2f} h\nz (m)", fontsize=11, fontweight="bold")
        axes[0, 0].legend(fontsize=6)
        fig.suptitle(title + " (rows: times, columns: points)")
        fig.tight_layout()
        path = os.path.join(OUT_DIR, filename)
        fig.savefig(path, dpi=110)
        print("wrote", path)

    cpuStyle = dict(color="k", ls="-", lw=1.2)
    gpuStyle = dict(color="r", ls="--", lw=1.0)
    profileFig([(uC, uG, cpuStyle, gpuStyle, "CPU u", "GPU u"),
                (vC, vG, dict(color="k", ls=":", lw=1.0), dict(color="tab:orange", ls=":", lw=1.0),
                 "CPU v", "GPU v")],
               "u, v (m/s)", RUN + " horizontal velocity profiles", PREFIX + "_velocity.png")
    profileFig([(wC, wG, cpuStyle, gpuStyle, "CPU", "GPU")],
               "w (m/s)", RUN + " vertical velocity profiles", PREFIX + "_wvel.png")
    if hasSalt:
        profileFig([(sC, sG, cpuStyle, gpuStyle, "CPU", "GPU")],
                   "salinity", RUN + " salinity profiles", PREFIX + "_salinity.png", minXMax=0.01)
        cpuS.close()
        gpuS.close()

    for ds in (cpu2d, gpu2d, cpuZ, gpuZ, cpuU, gpuU, cpuV, gpuV, cpuW, gpuW):
        ds.close()

    print("wrote", os.path.join(OUT_DIR, PREFIX + "_map.png"))
    print("wrote", os.path.join(OUT_DIR, PREFIX + "_eta.png"))
    plt.show()


if __name__ == "__main__":
    main()
