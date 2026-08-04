# CUDA Functions Module (`python/cuda_functions.py`)

The `cuda_functions.py` module serves as the primary Python interface for executing CUDA kernels defined in `kernels/meshdata.cu`. It handles the compilation, linking, and high-level management of CUDA kernels, abstracting away the complexities of `cupy.RawModule` and manual kernel launches.

## Overview

This module is responsible for:
1.  **Kernel Management**: Managing the lifecycle of the CUDA module via `CudaSourceManager`.
2.  **Kernel Loading**: Dynamically loading various kernel functions (e.g., `populateMesh`, `computeVerticalVelocity`) and binding them to Python-callable wrappers.
3.  **Kernel Dispatching**: Providing a standardized way to launch kernels based on their logical "type" (e.g., `global`, `nodes`, `elements`, `sides`).

## Core Functionality

### Initialization

The module uses two main initialization stages, typically called from `main.py`:

#### 1. `initializeCudaKernelManager(_project_root)`
This function initializes the `CudaSourceManager` and compiles the primary CUDA source file (`kernels/meshdata.cu`). 

- **Process**:
    - Sets the project root.
    - Uses `CudaSourceManager` to preprocess the `.cu` file (resolving `#include` directives).
    - Compiles the source into a `cupy.RawModule`.
    - Specifically extracts and wraps "global" type kernels: `createMesh`, `getStructSizes`, and `initTidalHarmonicConstituents`.
- **Global Kernel Concept**: A "global" kernel is specialized. Unlike other kernels that iterate over mesh components, a global kernel is launched with exactly **one thread**. Its purpose is to perform one-time setup, such as initializing global device variables or allocating memory.

#### 2. `initializeCudaFunctions(sizes)`
Once the module is loaded, this function iterates through a large list of available CUDA kernels and wraps them using the `CudaKernelWrapper`.

- **Parameters**:
    - `sizes` (dict): A dictionary containing the counts for different mesh entities (e.g., `{'nodes': 100, 'elements': 50, ...}`). These sizes are used to determine the grid dimensions for the kernel launch.
- **Kernel Types**: The wrapper uses the `kernel_type` to decide how many threads to launch:
    - **`global`**: 1 thread per block, 1 block total. Used for initialization.
    - **`nodes`**, **`elements`**, **`sides`**, **`nodeBnd`**, **`sideBnd`**: Standard parallel execution using `THREADS_PER_BLOCK` (256) and calculated number of blocks to cover all entities of that type.

## Kernel Function Library

The module exposes a vast array of functions. They can be categorized by their execution pattern:

### 1. Initialization & Setup (Global Type)
Kernels that run once with a single thread.
| Function | Purpose |
| :--- | :--- |
| `createMesh` | Initializes mesh structures on the GPU. |
| `getStructSizes` | Retrieves the sizes of various GPU structures. |
| `initTidalHarmonicConstituents` | Sets up tidal constituents. |
| `computeTimeRamps` | Computes time-dependent values. |
| `computeBedDeformation` | Computes bed deformation based on global parameters. |
| `initializeScalars` | Initializes scalar values on the device. |

### 2. Mesh Population (Entity-based Type)
Kernels that iterate over mesh components to transfer data from host to device or compute initial states.
| Function | Kernel Type | Description |
| :--- | :--- | :--- |
| `populateMesh` | `global` | (Actually defined as global in code, likely handles full mesh setup) |
| `populateNodes` | `nodes` | Populates node-specific data. |
| `populateSides` | `sides` | Populates side-specific data. |
| `populateElements` | `elements` | Populates element-specific data. |
| `populateBoundaryNodes` | `nodeBnd` | Populates boundary node data. |

### 3. Physics & Computation (Entity-based Type)
The core simulation kernels. These are heavily parallelized across nodes, elements, or sides.
| Function | Kernel Type | Description |
| :--- | :--- | :--- |
| `computeVerticalVelocity` | `elements` | Computes vertical velocity for each element. |
| `computeNodalVerticalVelocities1`| `elements` | Computes nodal vertical velocity from elements. |
| `computeHorizontalViscosity` | `sides` | Computes horizontal viscosity at sides. |
| `updateElevation` | `nodes` | Updates the water surface elevation at nodes. |
| `solveMomentumEquations` | `sides` | Solves momentum equations at sides. |
| `computeDensityGradient` | `elements` | Computes density gradients within elements. |
| `computeBottomDragVegetation` | `nodes` | Calculates drag forces due to vegetation at nodes. |

### 4. Data Access & Utilities
Functions to retrieve data or perform low-level memory operations.
| Function | Kernel Type | Description |
| :--- | :--- | :--- |
| `getElementVarAddr` | `global` | Retrieves address of a variable for an element. |
| `getNodeVarAddr` | `global` | Retrieves address of a variable for a node. |
| `getSideVarAddr` | `global` | Retrieves address of a variable for a side. |
| `printMatrix` | `global` | Prints a matrix for debugging. |
| `returnMatrix` | `global` | Returns a matrix from device to host. |

## Relationships

- **`main.py`**: The entry point that calls `initializeCudaKernelManager` and `initializeCudaFunctions`.
- **`iop/cudamodule.py` (`CudaSourceManager`)**: The engine that handles the actual compilation and pre-processing of the CUDA source.
- **`iop/cudaKernelWrapper.py` (`CudaKernelWrapper`)**: The wrapper class used by this module to encapsulate kernel launch logic and multi-GPU support.
- **`kernels/meshdata.cu`**: The source file containing the actual CUDA C++ implementations of all functions listed above.
- **`kernels/constants.h`**: A dynamically generated header file containing simulation constants (generated via `processConstants`) that is included in the CUDA compilation process.

## Implementation Notes

- **Multi-GPU Support**: The `CudaKernelWrapper` (invoked by these functions) is designed to launch the same kernel across all available GPUs in the system, provided the data is replicated or managed correctly on each device.
- **Memory Management**: Many kernels work with `cupy` arrays. The `sizes` dictionary passed during initialization is critical for defining the grid dimensions of the kernels.
