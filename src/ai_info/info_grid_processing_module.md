# Grid Processing Module (`read_and_process_grid.cpp`)

This module provides high-performance C++ implementations for reading and processing complex unstructured grids (horizontal and vertical). It is designed to be called from Python via **pybind11**, where it is compiled on the fly (or cached). 

The use of C++ is critical for this module because the I/O and geometric processing tasks (such as building adjacency lists and performing METIS partitioning) involve heavy computational loads and high-volume I/O that would be prohibitively slow in native Python.

## Overview

The module handles:
- Reading `.gr3` (horizontal grid) files.
- Reading `.ll` (longitude/latitude) files.
- Reading vertical grid files (`vgrid.in`) supporting both localized sigma and SZ coordinates.
- Reading wind rotation angle files.
- Computing connectivity structures: node-to-element, element-to-side, node-to-node, and element-to-element adjacency.
- Mesh partitioning using the **METIS** library.
- Reordering elements counter-clockwise around nodes.

## Key Functions

### `readHGrid`
Reads a horizontal grid file (`.gr3`) and associated longitude/latitude data.

**Signature:**
```cpp
std::tuple<
    py::array_t<int32_t>,    // i34: Number of nodes per element
    py::array_t<int32_t>,    // partElements: Partition index for each element
    py::array_t<int32_t>,    // partNodes: Partition index for each node
    py::array_t<int32_t>,    // partSides: Partition index for each side
    py::array_t<int32_t>,    // elnode: Global indices of nodes forming each element
    py::array_t<int32_t>,    // elside: Global indices of sides forming each element
    py::array_t<double>,     // nodeCoords: [x, y, z] coordinates of nodes
    py::array_t<double>,     // nodeCoords_ll: [lon, lat] coordinates of nodes
    py::array_t<int32_t>,    // ic3: Mapping from element side to neighbor element
    py::array_t<int32_t>,    // isdel: Mapping from side to adjacent elements
    py::array_t<int32_t>,    // isidenode: Mapping from side to its nodes
    py::array_t<int32_t>,    // isbnd: Boundary type for each node
    py::array_t<int32_t>,    // nond: Number of nodes per open boundary segment
    py::array_t<int32_t>,    // nlnd: Number of nodes per land boundary segment
    py::array_t<int32_t>,    // iond: Node indices for open boundary segments
    py::array_t<int32_t>,    // ilnd: Node indices for land boundary segments
    py::tuple,               // nodeNeighElems: ArrayOfRows (node -> elements)
    py::tuple                // nodeNeighNodes: ArrayOfRows (node -> nodes)
> readHGrid(const std::string& filename, const std::string& filename_ll, idx_t nPartitions, int ncor)
```

**Parameters:**
- `filename`: Path to the `.gr3` file.
- `filename_ll`: Path to the longitude/latitude file.
- `nPartitions`: Number of partitions to use for METIS.
- `ncor`: Flag to indicate if longitude/latitude data should be read.

---

### `readVGrid`
Reads vertical grid information, supporting two coordinate systems.

**Signature:**
```cpp
std::tuple<
    py::array_t<int32_t>,     // ivcor, nvrt (vertical coord type and number of levels)
    py::array_t<double>,       // ztot (depth at each level for SZ)
    py::array_t<int32_t>,      // kb (vertical level index for each node for localized sigma)
    py::array_t<double>,       // h_c, h_s, theta_b, theta_f (SZ parameters)
    py::array_t<double>,       // sigma_lcl (localized sigma coordinates)
    py::array_t<double>,       // sigma (sigma levels for SZ)
    py::array_t<double>,       // cs (transformed sigma)
    py::array_t<double>        // dcs (derivative of transformed sigma)
> readVGrid(const std::string& filename, int nNodes)
```

**Parameters:**
- `filename`: Path to the vertical grid file.
- `nNodes`: Number of nodes in the mesh.

---

### `readHGridData`
Reads scalar data associated with nodes from a `.gr3` file.

**Signature:**
```cpp
py::array_t<double> readHGridData(const std::string& filename, int nNodes, int nElements)
```

---

### `readWindRotationAngle`
Reads wind rotation angles for elements.

**Signature:**
```cpp
py::array_t<double> readWindRotationAngle(const std::string& filename, int nElements)
```

## Internal Helper Functions

### Grid Connectivity Construction
- `buildListsOfElementNeighbors`: Builds the mapping from nodes to the elements they belong to.
- `buildListsOfNodeNeighbors`: Builds the mapping from nodes to their neighboring nodes (nodes sharing an element).
- `buildAdjacencyList`: Creates an adjacency list for elements (element-to-element) for METIS partitioning.
- `buildElementSideElementTable`: Builds the `ic3` table which maps an element's side to its adjacent element.
- `reorderElementsCounterClockwise`: Reorders the elements around a node in a consistent CCW order.
- `reorderPartitionNodes`: Uses METIS to reorder nodes within a specific partition for locality.
- `partitionNodesAndSides`: Assigns nodes and sides to partitions based on the element partitions.

### Python Interop
- `toPyArray<T>` (overloads for `std::vector`, `Array2D`, `ArrayOfRows`): Converts C++ data structures into NumPy-compatible `py::array_t`.

## Dependencies
- **pybind11**: For Python bindings and NumPy interoperability.
- **METIS**: For graph partitioning and reordering.
- **Custom Headers**: `../include/misc.h`, `../include/grid_helper_types.h` for `Array2D`, `ArrayOfRows`, and utility functions like `progressBar`.
