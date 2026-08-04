### 1. Class Descriptions and Relationships
The mesh is represented as a graph of interconnected geometric and physical entities:
*   **`Node`**: Represents a discrete point in space. It stores nodal-centered variables (e.g., elevation `eta`, pressure) and maintains connectivity to the `Elements` and other `Nodes` that share it.
*   **`Element`**: Represents a computational cell (triangle or quadrilateral). It is the primary unit for integration and volume-averaging, holding pointers to its constituent `Nodes`, `Sides`, and neighboring `Elements`. It contains shape functions for interpolation.
*   **`Side`**: Represents the edge between two nodes and the interface between two elements. In a Finite Volume context, the `Side` is where fluxes (e.g., normal velocity `u2`) are computed. It links two `Nodes` and two `Elements`.
*   **`NodeBoundary` & `SideBoundary`**: These are "decorator" classes. If a `Node` or `Side` is on a domain boundary, it holds a pointer to one of these objects, which store boundary-specific parameters (e.g., time-history of elevation) and link back to a `BoundarySegment` defining the boundary type (e.g., "land" or "open").
*   **`Mesh`**: The top-level container. It manages memory allocation, global-to-local mapping, and the partitioning of the mesh across multiple GPUs.

**Relationship Summary:** 
`Mesh` $\rightarrow$ owns arrays of `Node`, `Element`, `Side`, `NodeBoundary`, and `SideBoundary`. `Element` $\leftrightarrow$ `Node` and `Element` $\leftrightarrow$ `Side` (Many-to-Many). `Side` $\leftrightarrow$ `Node` (Many-to-Many). `Node` $\rightarrow$ `NodeBoundary` and `Side` $\rightarrow$ `SideBoundary` (1:0 or 1).

### 2. Relation to the Numerical Model
This architecture directly maps to a **Finite Volume / Finite Element** numerical scheme:
*   **Nodal values** (`Node`) are used for state variables naturally defined at points.
*   **Face values** (`Side`) are used for fluxes; the simulation computes fluid movement across a side based on gradients between adjacent elements.
*   **Cell values** (`Element`) are used for volume-integrated quantities and conservation laws.
*   **Boundaries** allow for external forcings (tides, wind) to be applied only where necessary, without bloating the memory of internal entities.

### 3. CUDA and High-Performance Design
#### Pointers vs. Indices
The codebase uses direct pointers (`Node*`, `Element*`) instead of integer indices. In a CUDA environment, this is preferred because:
*   **Reduced Indirection**: Pointers avoid the calculation `BaseAddress + (Index * ElementSize)`, allowing the GPU to jump directly to the memory address.
*   **Faster Graph Traversal**: The solver frequently "walks" the mesh (e.g., `Element` $\rightarrow$ `Side` $\rightarrow$ `Neighbor Element`). Pointers turn this into a series of direct dereferences rather than repeated array lookups, reducing latency and instruction count.

#### Multi-GPU Systems
The `Mesh` class supports multi-GPU scalability via **Partitioning**:
*   **Pointer Arrays**: The use of double pointers (e.g., `Node** nodes`) allows the mesh to be split into partitions, where each pointer in the array points to the local memory of a specific GPU.
*   **Local Memory Management**: Each GPU manages its own `localNodeData`, `localElemData`, etc., while the `Mesh` structure provides the framework to access data across partitions (facilitated by Unified Memory or P2P transfers).
*   **Workload Balance**: The `partNum` identifies which slice of the global domain a GPU is processing, ensuring parallel execution across the cluster.

