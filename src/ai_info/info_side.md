The `Side` class (implemented as a `struct` in `kernels/meshdata.h`) is a fundamental component of the mesh topology, representing an edge shared between two elements (or a boundary edge). It serves as the primary interface for calculating fluxes, gradients, and performing boundary condition updates.

### **Overview**
In the simulation's unstructured grid, a `Side` acts as the connective tissue between `Element` and `Node` objects. While nodes store point-based data and elements store cell-centered data, the `Side` class manages data and computations associated with the faces of the elements, which is critical for the Finite Volume Method (FVM) used in the solver.

---

### **Data Members**

#### **Topology & Connectivity**
*   **`int idx`**: The global unique identifier for the side.
*   **`Element *elems[2]`**: Pointers to the two elements that share this side. For boundary sides, `elems[1]` is `nullptr`.
*   **`Node *n1, *n2`**: Pointers to the two nodes that define the endpoints of the side.
*   **`SideBoundary *sideBnd`**: A pointer to a `SideBoundary` object containing boundary-specific metadata. This is non-null only if the side is a boundary edge.

#### **Geometry**
*   **`double len`**: The physical length of the side.
*   **`Vector2 center`**: The 2D $(x, y)$ coordinates of the side's midpoint.
*   **`double zc[nvrt]`**: The vertical $z$-coordinates at the side center for each vertical level.
*   **`Vector2 n`**: The unit normal vector of the side.
*   **`int kb`**: The bottom index of the side.

#### **Fields**

*   **`Vector2 u2[nvrt]`**: Horizontal components of the velocity vector at every level.


#### **State**
*   **`bool isDry`**: A flag indicating if the side is currently "dry" (above the water surface), used for wetting and drying logic.

---

### **Key Methods**

#### **1. Topology Queries**
*   **`isBnd()`**: Returns `true` if the side is on the mesh boundary (`sideBnd != nullptr`).
*   **`isInternal()`**: Returns `true` if the side is shared by two elements (`elems[1] != nullptr`).
*   **`neigh(int index)`**: Returns a neighboring side based on the provided index, used for traversing the mesh.

#### **2. Geometric Interpolation**
The class provides parametric access to points along the edge using a parameter $t \in [0, 1]$:
*   **`x(t)`, `y(t)`, `z(t)`**: Linear interpolation between nodes `n1` and `n2` to find the coordinate at position $t$.
*   **`p2_t(t)` / `p3_t(t)`**: Helper methods that return the interpolated point as a `Vector2` or `Vector3`.
*   **`lon()`, `lat()`**: Returns the average longitude and latitude of the two endpoints.
*   **`centerZ()`**: Returns the average vertical coordinate of the two endpoints.

#### **3. Field Averaging**
`Side` implements several helper methods to compute the average value of a nodal field at the side center. These typically follow the pattern `0.5 * (n1->field + n2->field)`:
*   **Elevation/Bathymetry**: `eta1()`, `eta2()`, `depth()`.
*   **Physics/Friction**: `grav3()`, `Cd()`, `veg_h()`, `veg_nv()`.
*   **Dynamics**: `w2(int k)` (vertical velocity at level $k$).

#### **4. Advanced Numerical Operations**
*   **`grad2(...)`**: Computes the gradient of a specified field at the side center by averaging the gradients from the two adjacent elements.
*   **`tau()`**: Computes the stress tensor ($\tau$) at the side center, handling necessary projections between nodal and side frames.
*   **`btrack(...)`**: A complex backtracking method used for Lagrangian particle tracking or stability analysis, tracing a path backward in time from the side center.

---

### **Relationships to Other Classes**

| Class | Relationship | Description |
| :--- | :--- | :--- |
| **`Element`** | Many-to-Many | Each `Element` maintains an array of 3 or 4 `Side` pointers. Each `Side` knows its two adjacent `Element`s. |
| **`Node`** | Many-to-One | Each `Side` is defined by exactly two `Node` objects. |
| **`SideBoundary`** | One-to-One (Optional) | A `Side` may point to a `SideBoundary` object to handle boundary condition logic (e.g., open vs. closed boundaries). |
| **`Frame`** | Associated | Used via `Side` methods to transform vectors between global coordinates and local side-aligned frames. |

