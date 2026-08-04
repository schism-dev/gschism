# schism_step.py

`schism_step.py` is the main loop orchestrator for the SCHISM (Semi-implicit Cross-scale Hydroscience Integrated System Model) simulation. It coordinates the execution of various CUDA kernels and Python functions to perform a single time step of the simulation, covering processes such as tidal potential calculation, meteorological forcing, wave equation solving, momentum equation solving, and bed deformation updates.

## Functions

### `solveWaveEquation(solver, maxIter=1000, tol=1e-12)`

Solves the wave equation using the Conjugate Gradient method. It relies on an external solver object (likely a `PCGSolver` instance) and interfaces with the global simulation state.

**Parameters:**
| Name | Type | Description |
| :--- | :--- | :--- |
| `solver` | `PCGSolver` | An instance of the PCG solver used to solve the linear system. |
| `maxIter` | `int` | Maximum number of iterations for the solver. Defaults to 1000. |
| `tol` | `float` | Relative tolerance for convergence. Defaults to 1e-12. |

**Returns:**
| Type | Description |
| :--- | :--- |
| `int` | The number of iterations performed by the solver. |

**Notes:**
- Assumes all arrays are already stored in the `Node` structure as doubles.
- Uses global variables `N` and `NN` (though not explicitly passed, they are used for context).

---

### `computeMeteorologicalForcing(opt, time, sflux)`

Applies meteorological forcing (wind, pressure, etc.) to the simulation based on time-dependent data from an `sflux` dataset.

**Parameters:**
| Name | Type | Description |
| :--- | :--- | :--- |
| `opt` | `Object` | Configuration/option object containing parameters like `nws` (number of wind/surface forcings) and `start_year/month/day/hour`. |
| `time` | `float` | Current simulation time in seconds. |
| `sflux` | `Object` | An object (likely from `iop.readSflux`) that provides time-varying meteorological data. |

**Details:**
- Converts simulation time to Julian Day.
- Retrieves wind, pressure, and temperature data from `sflux`.
- Updates the GPU via `cufunc.setWinds`.

---

### `schism_step(it, namelist, bctides, sflux, solver)`

The core function that executes a single simulation time step. It calls a sequence of CUDA kernels and Python utilities to update the physical state of the model.

**Parameters:**
| Name | Type | Description |
| :--- | :--- | :--- |
| `it` | `int` | The current iteration (time step) index. |
| `namelist` | `Object` | A structure containing simulation configurations (`core` and `opt` subgroups). |
| `bctides` | `Object` | Object containing tidal boundary condition parameters (frequencies, amplitudes, etc.). |
| `sflux` | `Object` | Interface for reading surface flux data. |
| `solver` | `PCGSolver` | The iterative solver used for the wave equation. |

**Execution Flow:**
1. **Initialization & Scaling:** Scales gravity based on tidal loading coefficients and computes time ramps for various forcing types.
2. **Forcing Application:** 
   - Computes bed deformation.
   - Computes Earth tidal potential and loading tides.
   - Applies meteorological forcing (wind/pressure).
   - Sets boundary conditions.
3. **Hydrodynamics (Momentum & Mass):**
   - Computes bottom drag and vegetation drag.
   - Applies horizontal viscosity (Shapiro filter if enabled).
   - Uses `backtrackAllNodes` (upwind/backtracking) for momentum.
   - Computes vertically integrated velocity.
   - Computes density gradients and baroclinic terms.
4. **Wave Equation Solver:**
   - Sets up the coefficient matrix and RHS.
   - Performs the roll-up of elevation.
   - Solves the wave equation using the provided `solver`.
   - Updates elevation.
5. **Momentum Updates:**
   - Computes horizontal gradients.
   - Solves momentum equations and applies boundary conditions.
6. **Sub-grid / Wet-Dry Dynamics:**
   - Computes vertical velocity.
   - Updates Z-coordinates (bed deformation).
   - Manages dry/wet element and node transitions (`ul0` functions).
   - Updates sides and nodal velocities.

**Dependencies:**
- `python.cuda_functions` (as `cufunc`): The primary interface for GPU-accelerated computations.
- `python.pcgsolver` (as `PCGSolver`): The solver engine.
- `iop.readSflux` and `iop.writeOutput`: I/O operations.
- `cupy` and `numpy`: Array manipulations.
- `debug.monitor`: For debugging and monitoring.
