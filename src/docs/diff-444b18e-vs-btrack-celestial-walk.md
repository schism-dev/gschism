# Differences: `444b18e` vs `btrack-using-celestial-walk`

**Baseline (old):** `444b18e` — *Works for ivcor==1…*  
**Branch (new):** `btrack-using-celestial-walk`  
**Generated:** 2026-05-29  

This document summarizes changes relevant to **GPU vs CPU / old vs new** comparisons, especially for **`ics==1`, `ibc==1`, `iveg==0`**.

---

## How to regenerate artifacts

From the repo root (`gschism/`):

```bash
# Human-readable summary (this file)
# Edit src/docs/diff-444b18e-vs-btrack-celestial-walk.md

# Full unified diff (core sources, ~240 KB)
git diff 444b18e..btrack-using-celestial-walk -- \
  src/kernels/meshdata.cu src/kernels/meshdata.h \
  src/python/schism_step.py src/kernels/constants.h \
  src/iop/writeOutput.py \
  > src/docs/diff-444b18e-vs-btrack-celestial-walk.patch

# Full tree stat
git diff 444b18e..btrack-using-celestial-walk --stat

# Full patch (everything, very large)
git diff 444b18e..btrack-using-celestial-walk > src/docs/full.patch

# PDF/HTML: open the .patch in VS Code / Cursor, or use:
#   git difftool 444b18e..btrack-using-celestial-walk
```

**Files in this folder:**

| File | Contents |
|------|----------|
| `diff-444b18e-vs-btrack-celestial-walk.md` | This summary |
| `diff-444b18e-vs-btrack-celestial-walk.patch` | Unified diff (kernels + driver + output) |

---

## Commits on the branch (since `444b18e`)

```
70a4860, 6452eb9, 72fad4c, b9d71c4, 44b1330, add9ffd, ecf3f5c, cc3a0e5,
3a8a6de, e5ac1fe, 4763797, 00f2e03, 8457f03, bff9ee9, 4455ea9, 5f844fd,
2d9b4aa, 5a5b3ca, fa31f73, 7f6e363, 5f551f2, …
```

Notable messages: celestial walk / btrack, **ETP & vertical velocity fixes** (`e5ac1fe`), **ics==2**, element depth fix (`5f844fd`).

---

## Files changed (overview)

| Area | Files | Notes |
|------|-------|--------|
| Kernels | `meshdata.cu`, `meshdata.h` | Largest physics diff |
| Timestep | `schism_step.py` | Call order, sflux time, BC areas |
| I/O | `writeOutput.py` | Output / plotting hooks |
| Debug | `debug/monitor.py`, `newmonitor.py`, `render.cpp` | Visualization |
| Build | `constants.h`, `cuda_functions.py`, `cudamodule.py` | Compile-time / bindings |

---

## Not changed (for typical namelist)

With **`ibc==1`, `iveg==0`, `ics==1`, `ihorcon==0`, `itur==0`**:

| Topic | Status |
|-------|--------|
| Vertical viscosity (`dfv`) | Still constant `dfv0` in `populateNodes`; no `computeTurbulence` |
| Horizontal viscosity (`nuGradU`) | Still zeroed when `ihorcon==0` |
| Baroclinic density (`bcc`, `gradRho` integration) | Off (`ibc==1`) |
| Hwang wind stress (`iwind_form==1`) | Same formulas |
| `computeNodalHorizontalVelocities` (ics==1) | Same `u = s2+s3-s1` |
| `solveMomentumEquations` 3D core | Same tridiagonal structure |

---

## High-impact differences (velocity / wetting)

### 1. Backtracking (`Side::btrack` → `celestialWalk`)

- Old: element-to-element walk + `trackThroughElementsAndLevels`.
- New: **Leaded Celestial Walk** (`celestialWalk`), boundary intersection logic, `sincosf(-wind_rotate_angle)` for wind.
- Affects **`btu`** → momentum RHS; largest intentional change if btrack is not held fixed.

### 2. `computeVerticalVelocity` — **strong effect on `w` / vertical profiles**

| Item | Old (`444b18e`) | New (branch) |
|------|-----------------|--------------|
| `elem->w[kb]` | `(av_bdef2-av_bdef1)/dt - dot(gradh,u2meanBottom)` | `bedVelocity + bottomW`, `bedVelocity=0` for `imm==0` |
| `bottomFlux` | `-dot(gradh, u2meanBottom)` | **`0`** (no bed-face flux) |
| Side `u` in flux | `u2[l]`, `u2[l+1]` | `u2[max(l,kbs)]`, `u2[max(l+1,kbs)]` |
| Side area | `side->area(l)` | **`side->areaClamped(l)`** |
| Top `w` | combined flux / `(nz*area)` | split horizontal vs vertical flux |

→ Changes **`elem->w`** and **`node->w2`** via `computeNodalVerticalVelocities*`.  
→ **Horizontal `u2` unchanged** if `side->u2` unchanged.

### 3. `computeElemAverageGhat` — ETP/pressure finish-off loop

- Fortran: `do j=1,i34` on **all** sides.
- Old CUDA: `for (j=1; j<i34; j++)` — **skipped side index 0** (porting bug).
- New CUDA: `for (j=0; j<i34; j++)` — matches Fortran.
- Active for **`ibc==1`** (tidal potential + `pr`); **not** baroclinic density.
- Can change **`ghat1`** → elevation → `htot`, `chi`, implicit momentum.

### 4. `interpolateVelocities` (btrack / element interp)

- Old: `L = kb` for vertical bracketing search.
- New: `L = kb - 1` (comment: element vs node bottom level).
- New: `getNodalValues(level, …)` uses **`max(level, node->kb)`** for nodal fields.

### 5. `populateBoundaryNodes` (tides)

- Tide index fix for Fortran-order `emo(nope, mnond, nbfr)` (CuPy `order='F'`).
- Old: `idx2 = iseg*mnond*nbfr + i + jfr*mnond` (wrong).
- New: `idx2 = iseg + i*nope + jfr*nope*mnond`.

### 6. `schism_step.py` (driver)

- `julian_day += opt.utc_start/24.0` in meteorological forcing.
- `resetBoundaryAreas()` / `computeBoundaryAreas()` before flow BCs.
- Debug hooks (`newmonitor.show()` at large `it`, etc.).
- `iupwind_mom != 0` guard.

### 7. `ics==2` only (not your case)

- Frame projections in `ghat1`, Shapiro, nodal interp, `computeBaroclinicTerm`, etc.
- Many `assert(false)` stubs **implemented** on new branch.

---

## Baroclinic / viscosity (your parameters)

| Piece | `ibc==1` behavior | Branch delta |
|-------|-------------------|--------------|
| `computeDensityGradient` | Early return | + `assert(ibc==1)` only |
| `computeBaroclinicTerm` | `if (ibc==1) return` | `ics==2` projection only |
| `ghat1` density block | Compiled out | Loop bound `k>=kb+1` vs `k>kb+1` **inactive** |
| `dfv` / turbulence | Constant `dfv0` | None |

---

## Wind stress (`nws==2`, `ihconsv==0`, `iwind_form==1`)

- Hwang formula: **same** on both versions.
- CUDA rotates with `cRot`/`sRot`; Fortran uses `windx/windy` — frame must match inputs.

---

## Known issues (both versions unless noted)

| Issue | Where | Notes |
|-------|--------|------|
| `kin = k - kb - 1` at `k==kb` | `solveMomentumEquations` | `soln[-1]`; often masked by `u2[kb]=u2[kb+1]` if `Cd==0` |
| Nodal interp without `max(k,side->kb)` | `computeNodalHorizontalVelocities1` | Both versions; near `kb` transitions |
| `iter>39` break in btrack | New branch (debug) | Truncates backtracking if present |

---

## Suggested A/B tests (velocity profiles)

1. Compare **`side->u2[k]`** old vs new (isolates momentum + shapiro).
2. If (1) matches, compare **`node->u2[k]`** (nodal interp).
3. Compare **`elem->w` / `node->w2`** (isolates `computeVerticalVelocity`).
4. Hold btrack fixed: same `btu` input or disable celestial-walk path.
5. Compare elevation / `kb` if profiles differ only near wet/dry.

---

## P&P wind (`iwind_form==-1`) — not your case

Fortran `dragcmin`/`dragcmax` at 6 and 50 m/s ≡ CUDA `clamp(speed,6,50)` before P&P formula; equivalent for `Cd`.

---

## References

- Fortran `ghat1` loop: `do j=1,i34` for both node/side pass and ETP finish-off (`do j=1,i34 !side`).
- Conversation context: CPU/GPU compare for `ics==1`, `ivcor==1`, `kb` transitions.
