In the codebase, the geometric types for points, vectors, and frames are defined in `kernels/meshdata.h` and leverage CUDA's built-in vector types for efficiency.

### 1. Points and Vectors
Points and vectors are implemented as type aliases for CUDA's `double2` and `double3` structures. The distinction between a `Point` and a `Vector` is semantic rather than functional at the compiler level.

*   **2D Geometry**: Both `Point2` and `Vector2` are aliases for **`double2`**.
*   **3D Geometry**: Both `Point3` and `Vector3` are aliases for **`double3`**.

The codebase provides a comprehensive set of overloaded operators and helper functions for these types, including:
*   **Arithmetic**: Addition (`+`, `+=`), subtraction (`-`, `-=`), and scalar multiplication/division.
*   **Linear Algebra**: `dot` product, `cross` product (for 2D and 3D), `length`, `norm2`, and `normalize`.
*   **Utilities**: `rot90` (for 2D vectors) and conversion functions between 2D and 3D representations (e.g., `toVector2`, `toPoint3`).

### 2. Frames
A **`Frame`** is a custom `struct` used to handle coordinate transformations between global Cartesian coordinates and local frames (typically East, North, Up).

*   **Representation**: To optimize memory, a `Frame` does not store the full $3 \times 3$ rotation matrix. Instead, it stores 5 coefficients:
    *   **Column 0 (East)**: `t00`, `t10` (the $z$-component `t20` is implicitly $0$).
    *   **Column 1 (North)**: `t01`, `t11`, `t21` (a normalized vector on the oblate ellipsoid).
    *   **Column 2 (Up)**: Not stored; it is recomputed on demand as the cross product of the East and North columns via the `getCol2()` method.

*   **Key Functionalities**:
    *   **`update(lat, lon)`**: Recalculates the frame's basis vectors based on provided latitude and longitude.
    *   **`G2L(Vector3)`**: Transforms a vector from the **Global** frame to the **Local** frame.
    *   **`L2G(Vector3)`**: Transforms a vector from the **Local** frame to the **Global** frame.
    *   **`projHorVec(Vector2, Frame)`**: Transforms a horizontal vector from the current frame to another frame and projects it onto the new local horizontal plane.

