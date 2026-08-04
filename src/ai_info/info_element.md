The `Element` class (defined as a `struct` in `kernels/meshdata.h`) is a central component of the mesh representation. It represents a computational cell—either a triangle or a quadrilateral—and manages both the geometric properties of the cell and the physical state variables distributed across its vertical levels.

### Overview
The `Element` class serves as the primary unit for finite volume/element calculations. It maintains connectivity to its constituent nodes and sides, tracks neighboring elements, and stores vertically-indexed physical data (such as density and tracers). It also provides the mathematical machinery for shape function interpolation and coordinate transformations between global and local frames.

---

### Data Members

#### 1. Geometry and Topology
*   **Connectivity**:
    *   `int idx`: Unique identifier for the element.
    *   `int i34`: Number of sides/nodes (3 for triangles, 4 for quadrilaterals).
    *   `Side *sides[4]`: Pointers to the boundary sides of the element.
    *   `Node *nodes[4]`: Pointers to the vertices of the element.
    *   `Element *nghbElems[4]`: Pointers to adjacent elements.
*   **Spatial Properties**:
    *   `double area`: The horizontal area of the element.
    *   `double radius`: The equivalent radius of the element.
    *   `Point2 center`: The 2D coordinates of the element center.
    *   `double ssign[4]`: Signs associated with each side, used in flux calculations.
*   **Vertical Structure**:
    *   `int kb`: The bottom index of the element in the vertical coordinate system.
    *   `double depth`: The minimum depth across all nodes of the element.
    *   `double z[nvrt]`: The vertical coordinates (local frame, vertical up) for each level.
*   **Local Frame & Gradients**:
    *   `Frame frame`: (Conditional on `USE_FRAMES`) A local coordinate frame used to handle rotations and projections.
    *   `Vector2 gradPhi[4]`: Derivatives of the shape functions.
    *   `Vector2 gradEta`, `gradPr`, `gradETP`: Gradients of free surface elevation, pressure, and Earth Tidal Potential, respectively.

#### 2. Physical State
Most physical variables are stored as arrays of size `nvrt` (number of vertical levels), representing values at prism centers.
*   **Dynamics**:
    *   `double w[nvrt]`: Vertical velocity at element centers across levels.
    *   `Vector2 ghat1`: Element-averaged $\breve{G}$ (used in momentum equations).
    *   `double vsource`: Volumetric source/sink term $[m^3/s]$.
*   **Density**:
    *   `double rho[nvrt]`: Density at each prism center.
    *   `double rho_mean[nvrt]`: Mean density at each prism center.
    *   `CubicSpline rho_spl`: Cubic spline interpolation for density.
    *   `double rho_z[nvrt]`, `rho_zz[nvrt]`: First and second derivatives of density.
    *   `Vector2 gradRho[nvrt]`: Gradient of de-meaned density.
*   **Tracers**:
    *   `Tracers tracers[nvrt]`: Arrays of tracer values at each vertical level.
*   **State Flags**:
    *   `bool isDry`: Indicates if the element is currently dry.

---

### Member Functions

#### 1. Topology and Classification
*   `isTri()` / `isQuad()`: Returns whether the element is a triangle or quadrilateral.
*   `is2D()` / `is3D()`: Determines the dimensionality of the element based on the bottom index `kb`.
*   `nodeLocalIdx(Node *node)` / `sideLocalIdx(Side *side)`: Returns the local index (0–3) of a given node or side within the element.
*   `nxq(int index)`: A utility for wrapping indices based on the number of sides `i34`.

#### 2. Vertical Coordinate Management
*   `findLevel(double z1, int k)`: Searches for the vertical level index containing the coordinate `z1`, starting the search from level `k`.
*   `z_prism(int k)`: Returns the midpoint of vertical level `k`.
*   `Dz(int k)`: Returns the thickness of vertical level `k`.
*   `centerZ()`: Returns the vertical center coordinate of the element.

#### 3. Mapping and Interpolation
*   **Coordinate Transformation**:
    *   `outputCoords()`: Returns the coordinates for output purposes.
    *   `getLocalHorizontalNodePosition(int idxNode)`: Computes the local horizontal position of a node, ensuring the element is treated as horizontal.
*   **Shape Functions (for Quadrilaterals)**:
    *   `Phi(int idxNode, Point2 p)`: Evaluates the shape function $\Phi$ for a node at local point $p(\eta, \xi)$.
    *   `dPhi_dxi(...)` / `dPhi_deta(...)`: Evaluates the derivatives of the shape function.
    *   `GaussPoint(int idxNode)`: Returns the $(\eta, \xi)$ coordinates of the Gauss-Legendre quadrature points.
    *   `jacobianDet(Point2 p)`: Calculates the determinant of the Jacobian at local point $p$.

#### 4. State and Physics Utilities
*   `isDry2()`: A more up-to-date check for whether the element is dry, typically using current elevation.
*   `demeanedDensity(int k)`: Returns the difference between the actual density and the mean density at level `k`.
*   `isEdgeBoundary(Side *side)`: Checks if a specific side of the element is on the mesh boundary.
*   **Aggregates** (implemented in `meshdata.cu`):
    *   `interpolate(...)`: Interpolates values within the element.
    *   `mean(...)`: Calculates the average value across the element.
    *   `getNodalValues(...)` / `getSideValues(...)`: Retrieves values associated with nodes or sides.
