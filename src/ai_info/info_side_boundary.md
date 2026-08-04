The `SideBoundary` class (implemented as a `struct` in `kernels/meshdata.h`) is designed to store and manage boundary condition data specifically associated with the edges (sides) of the mesh. It acts as a bridge between a physical `Side` and its corresponding `BoundarySegment`, providing utility methods to interpolate boundary values from nodes to the side center.

### **Class Overview**
`SideBoundary` is used primarily in CUDA kernels to handle the enforcement of velocity and elevation boundary conditions at the edges of the computational domain. While `NodeBoundary` stores data at the vertices, `SideBoundary` facilitates the retrieval of these values at the midpoints of the sides.

### **Data Members**
- **`Side *side`**: A pointer to the `Side` object that this boundary information belongs to.
- **`BoundarySegment *segment`**: A pointer to the `BoundarySegment` defining the properties of this boundary region.
- **`Vector2 ubc[nvrt]`**: An array containing velocity boundary conditions for each vertical level (`nvrt`).

### **Key Methods**

#### **Interpolation and Weighting**
Because boundary values are often defined at nodes, `SideBoundary` provides a mechanism to calculate a weighted average at the side center.

- **`getBoundaryArrayValueWeight()`**
    - **Purpose**: Determines the interpolation weight $w$ for the first node (`n1`) of the side.
    - **Logic**: 
        - If both nodes belong to the same `BoundarySegment`, it returns `0.5` (simple average).
        - If only `n1` matches the segment, it returns `1.0`.
        - If only `n2` matches the segment, it returns `0.0`.
    - **Return**: A `double` representing the weight of `n1`. The weight of `n2` is implicitly $(1 - w)$.

#### **Value Retrieval (Interpolated)**
The following methods use the weighting logic above to compute values at the side center:

- **`etaMeanBnd()`**
    - Returns the arithmetic mean of the `eta_mean` values from the two connected nodes' `NodeBoundary` objects.
- **`umo(int jfr)`**
    - Returns the interpolated velocity $\mathbf{u}_{mo}$ at the side center for a given frame `jfr`.
- **`ufa(int jfr)`**
    - Returns the interpolated velocity $\mathbf{u}_{fa}$ at the side center for a given frame `jfr`.
- **`uth(int k)`**
    - Returns the interpolated velocity $\mathbf{u}_{th}$ at the side center for vertical level `k`.
- **`nuth(int k)`**
    - Returns the **normal component** (the $x$-component in the local boundary frame) of the interpolated velocity $\mathbf{u}_{th}$ at vertical level `k`.

### **Integration with Mesh Topology**
- **`Side` $\rightarrow$ `SideBoundary`**: Each `Side` struct contains a pointer `SideBoundary *sideBnd`. If the side is not a boundary, this pointer is `nullptr`.
- **`Node` $\rightarrow$ `SideBoundary`**: `NodeBoundary` objects can also maintain references to the `SideBoundary` objects associated with the sides connected to that node.
- **Initialization**: The `Mesh::populateSides` method is responsible for identifying which sides are boundaries (via `findSideBoundaryType`) and allocating/assigning `SideBoundary` objects from the `localSideBoundaryMem` pool.

