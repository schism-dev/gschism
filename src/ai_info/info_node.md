The `Node` class is a fundamental structural component of the simulation's mesh, representing a discrete point in 3D space. It serves as the primary storage for nodal variables (such as elevation and velocity) and maintains the topology required for spatial discretization and the solution of the governing equations.

### Class Overview
The `Node` class is implemented as a C++ `struct` designed for high-performance access within CUDA kernels. It stores both static geometric data and dynamic physical state variables.

---

### Data Members

#### 1. Geometry and Topology
These members define the node's position and its connectivity within the mesh.
- **`idx`**: Unique identifier for the node.
- **`p`**: `Point2` horizontal coordinates. If `ics == 2`, these are the X and Y components of the global 3D Cartesian coordinates.
- **`gz`**: Global vertical coordinate.
- **`lon`, `lat`**: Geographic longitude and latitude in radians.
- **`depth`**: Bathymetric depth at the node.
- **`z[nvrt]`**: Local vertical coordinates at each vertical level.
- **`kb`**: The bottom index of the node, used to identify the lowest active vertical level.
- **`nodes` / `numNodes`**: A list of pointers to neighboring nodes and the count of those neighbors.
- **`elements` / `numElements`**: A list of pointers to elements that share this node and the count.
- **`numSides`**: The number of sides connected to this node.
- **`frame`**: (Conditional: `USE_FRAMES`) A `Frame` object providing a local coordinate system for rotation and projection.

#### 2. Physical State
Dynamic variables updated during the simulation time-stepping.
- **`eta1`, `eta2`**: Free surface elevation at the previous (`eta1`) and current (`eta2`) timesteps.
- **`u2[nvrt]`, `w2[nvrt]`**: Horizontal (`u2`) and vertical (`w2`) velocity components at each vertical level.
- **`utmp[nvrt]`**: Temporary storage for velocity calculations.
- **`rho[nvrt]`**: Density at each vertical level.
- **`dfv[nvrt]`**: Vertical viscosity at each vertical level.
- **`tracers[nvrt]`**: Tracer concentrations at each vertical level.
- **`dav`**: Depth-averaged velocity vector.
- **`ETP`, `pr`, `grav2`**: Earth Tidal Potential, air pressure, and local gravity.

#### 3. Vegetation and Bottom Drag
Parameters used to calculate friction and resistance.
- **Vegetation**: Includes `veg_h` (canopy height), `veg_nv` (vegetation density), and `veg_alpha3D[nvrt]` (vertical distribution of vegetation).
- **Bottom Drag**: 
    - `Cd`: Bottom drag coefficient.
    - **Roughness Union**: A mutually exclusive set of parameters (`rmanning`, `rough`, or `drag`) used to define the bottom boundary layer.
    - `tau_bot`: The resulting bottom stress vector.

#### 4. Solver State
Variables used specifically by the linear system solver (e.g., PCG).
- **Matrix Coefficients**: `coeffDiag` (diagonal element) and `coeffMatrix` (pointer to coefficients for neighbor nodes).
- **RHS**: `qel` represents the right-hand side of the elevation equation.
- **PCG Variables**: `mP` (search direction), `mR` (residual), `mAp` (matrix-vector product), and `mX` (solution).

#### 5. Meteorological Forcing
Inputs from atmospheric models.
- **`u_air`, `p_air`, `t_air`, `q_air`**: Air velocity, pressure, temperature, and humidity.
- **`tau`**: The surface stress vector exerted by the wind.
- **`cRot`, `sRot`**: Cosine and sine of the mesh rotation angle.

---

### Methods

#### Geometric Utilities
- **`outputCoords()`**: Returns the node's horizontal position. It returns Cartesian `p` if `ics == 1`, or longitude/latitude (converted to degrees) if `ics == 2`.
- **`Dz(int k)`**: Calculates the vertical distance between level $k$ and $k-1$.
- **`htot2()`**: Returns the total water column height (`eta2 + depth`).
- **`zTopCanopy()`**: Computes the absolute vertical elevation of the top of the vegetation canopy.
- **`isDry2()`**: Evaluates if the node is "dry" based on the updated total water height compared to a threshold `h0`.

#### Physical Updates
- **`setWind(Vector2 &u)`**: Rotates the incoming wind vector `u` into the local node frame and calculates the resulting surface stress `tau` using parameterized drag coefficients (based on `iwind_form`).
- **`setFrame(double lat, double lon)`**: Updates the internal `Frame` object based on the provided geographic coordinates.
- **`computeBottomDragCd()`**: Updates the bottom drag coefficient `Cd` based on the chosen roughness model.

#### Solver Support
- **`clearMatrixRow()`**: Resets the diagonal and neighbor coefficients in the sparse matrix to zero.
- **`addToSparseMatrix(Node const *node, double value)`**: Adds a value to the matrix entry for the current node and the specified neighbor node.
- **`updateZCoordData()`**: Synchronizes and updates the vertical coordinate data.

---

### Relationships
- **`NodeBoundary`**: If a node is on the boundary, `nodeBnd` points to a `NodeBoundary` object containing boundary-specific conditions (e.g., tidal amplitudes, phases, and essential BC flags).
- **`Element` & `Side`**: Nodes are the vertices for `Element` and `Side` objects. Elements use nodal values for shape-function interpolation.

