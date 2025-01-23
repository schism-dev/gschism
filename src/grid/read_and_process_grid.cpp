#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <metis.h>
#include <cassert>  

#include "../include/misc.h"
#include "../include/grid_helper_types.h"

namespace py = pybind11; 

// 

// elnode[4, nElements] -> Global index to the nodes that form the element
// elside[4, nElements] -> Global index to the sides of the element
// i34[nElements] -> Number of nodes per element
// ic3[4, nElements] -> Global index to the element that shares the side 
// isdel[2, nSides] -> Global index to the elements that share the side

void buildVertexElements(const Array2D<int>& elnode, const std::vector<int>& i34, int nElements, int nNodes, 
                         ArrayOfRows<int>& vertexElements) {

    // Initialize arrays
    std::vector<size_t> nne(nNodes, 0);  // vertexConnectionCounts.

    // Count elements connected to each vertex
    for (int i = 0; i < nElements; i++) {
        int nVerts = i34[i];
        
        // Add one to connection count for each vertex in this element
        for (int j = 0; j < nVerts; j++) {
            int vertex = elnode(j, i);   
            nne[vertex]++;
        }
    }

    // Check for hanging nodes
    int nHanging = 0;
    for (int ip = 0; ip < nNodes; ip++) {
        if (nne[ip] == 0) {
            printFormatted("Hanging node found: %d", ip + 1);
            nHanging++;
        }
    }

    if (nHanging > 0) {
        char msg[256];
        snprintf(msg, sizeof(msg), "%d hanging nodes found", nHanging);
        throw std::runtime_error(msg);
    }

    // Create the vertex-to-elements mapping "array"
    vertexElements.rebuild(nne);

    // Build vertex-to-elements mapping
    std::fill(nne.begin(), nne.end(), 0);  // Reset connection counts
    for (int elemIdx = 0; elemIdx < nElements; elemIdx++) {
        progressBar("Building vertex-to-element mapping", elemIdx, nElements, ansi::YELLOW, 20);
        int nVerts = i34[elemIdx];

        // For each vertex in this element
        for (int i = 0; i < nVerts; i++) {
            int vertex = elnode(i, elemIdx);   

            // Add element to vertex's connection list
            const int j = nne[vertex];

            // Add element to vertex's connection list
            vertexElements(vertex, j) = elemIdx;

            nne[vertex]++;
        }
    }

    print("");
}


AdjacencyList buildAdjacencyList(const Array2D<int>& elnode, const std::vector<int>& i34, 
                                 const ArrayOfRows<int>& vertexElements, int nElements) {
    // First pass: count edges
    int ntedge = 0;  // total number of edges
    int mxnedge = 0;  // max number of local edges
    std::vector<int> tempAdjncy(1000, 0);  // temporary array for single element

    for (int ie = 0; ie < nElements; ie++) {
        progressBar("Building adjacency list", ie, nElements, ansi::MAGENTA, 20);
        int nedge = 0;  // number of local edges
        std::fill(tempAdjncy.begin(), tempAdjncy.end(), 0);  // reset array

        // For each node in element
        for (int j = 0; j < i34[ie]; j++) {
            int ip = elnode(j, ie);

            // For each element connected to this node
            for (size_t k = 0; k < vertexElements.rowSize(ip); k++) {

                const int je = vertexElements(ip, k);

                // Check if je is not ie and not already in tempAdjncy
                if (je != ie) {
                    bool found = false;
                    for (int m = 0; m < nedge; m++) {
                        if (tempAdjncy[m] == je) {
                            found = true;
                            break;
                        }
                    }
                    
                    if (!found) {
                        if (nedge >= 1000) {
                            throw std::runtime_error("Number of edges exceeds temporary array size (1000)");
                        }
                        tempAdjncy[nedge] = je;
                        nedge++;
                    }
                }
            }
        }
        ntedge += nedge;
        mxnedge = std::max(mxnedge, nedge);
    }
    print("");

    // Allocate arrays for graph
    AdjacencyList result;
    result.xadj.resize(nElements + 1, 0);
    result.adjncy.resize(ntedge, 0);
    result.adjwgt.resize(ntedge, 0);
    result.ntedge = ntedge;

    // Second pass: build edge list and weights
    result.xadj[0] = 0;

    for (int ie = 0; ie < nElements; ie++) {
        progressBar("Second pass: edge list and weights", ie, nElements, ansi::MAGENTA, 20);
        int nedge = 0;

        for (int j = 0; j < i34[ie]; j++) {
            int ip = elnode(j, ie);

            for (size_t k = 0; k < vertexElements.rowSize(ip); k++) {

                const int je = vertexElements(ip, k);

                if (je != ie) {
                    // Check if je is not already in adjncy
                    bool found = false;
                    for (int l = result.xadj[ie]; l < result.xadj[ie] + nedge; l++) {
                        if (result.adjncy[l] == je) {
                            result.adjwgt[l] = (2 * i34[je] - 3);
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        result.adjncy[result.xadj[ie] + nedge] = je;
                        result.adjwgt[result.xadj[ie] + nedge] = (2 * i34[je] - 1);
                        nedge++;
                    }
                }
            }
        }
        result.xadj[ie + 1] = result.xadj[ie] + nedge;
    }
    print("");

    return result;
}


void buildElementSideElementTable(const Array2D<int>& elnode, const std::vector<int>& i34,
                                  const ArrayOfRows<int>& vertexElements, 
                                  int nElements, int nNodes, 
                                  Array2D<int>& elside, Array2D<int>& isdel, Array2D<int>& ic3,
                                  int &nSides) {

    // Builds the table (ic3)that maps each side of an element to the element that shares the side


    for (int ie = 0; ie < nElements; ie++) {
        progressBar("Building element-side-element table", ie, nElements, ansi::CYAN, 20);
        for (int k = 0; k < i34[ie]; k++) {
            // Get nodes for this side
            const int loc_n1 = nxQ(0, k, i34[ie]);
            const int loc_n2 = nxQ(1, k, i34[ie]);
            const int n1 = elnode(loc_n1, ie);
            const int n2 = elnode(loc_n2, ie);

            // Look for matching element
            for (size_t l = 0; l < vertexElements.rowSize(n1); l++) {

                const int je = vertexElements(n1, l);

                // Check if this element shares n2
                if (je != ie) {
                    // Check if n2 exists in element je
                    bool found = false;
                    for (int m = 0; m < i34[je]; m++) {
                        if (elnode(m, je) == n2) {
                            found = true;
                            break;
                        }
                    }
                    
                    if (found) {
                        ic3(k, ie) = je;
                        break;
                    }
                }
            }

            // Check orientation
            int je = ic3(k, ie);
            if (je != -1) {
                for (int l = 0; l < i34[je]; l++) {
                    int n1_je = nxQ(0, l, i34[je]);
                    int n2_je = nxQ(1, l, i34[je]);

                    if (elnode(n1_je, je) == n1 && elnode(n2_je, je) == n2) {
                        throw std::runtime_error("Elem " + std::to_string(ie+1) + 
                                              " and " + std::to_string(je+1) + 
                                              " have opposite orientation");
                    }
                }
            }
        }
    }

    // Count global number of sides and build global element-side index table
    int ns = 0;

    printFormatted(ansi::blue("Building global element-side index table"));
    for (int ie = 0; ie < nElements; ie++) {
        for (int j = 0; j < i34[ie]; j++) {  // Visit each side associated with element ie
            if (ic3(j, ie) == 0 || ie < ic3(j, ie)) {  // New global side
                ns++;
                elside(j, ie) = ns;
                
                if (ic3(j, ie) != 0) {  // Old internal side
                    int je = ic3(j, ie);
                    int l = -1;
                    
                    // Find matching side in neighboring element
                    for (int k = 0; k < i34[je]; k++) {
                        if (ic3(k, je) == ie) {
                            l = k;
                            break;
                        }
                    }
                    
                    if (l == -1) {
                        throw std::runtime_error("AQUIRE_HGRID: Wrong ball info for elements " + 
                                              std::to_string(ie) + ", " + 
                                              std::to_string(j) + ", " + 
                                              std::to_string(ns));
                    }
                    elside(l, je) = ns;
                    // isdel(1, ns) = je;
                }
            }
        }
    }

    if (ns < nElements || ns < nNodes) {
        throw std::runtime_error("AQUIRE_HGRID: weird grid with ns < ne or ns < np: " + 
                               std::to_string(nNodes) + ", " + 
                               std::to_string(nElements) + ", " + 
                               std::to_string(ns));
    }


    // Build isdel
    nSides = ns + 1;
    isdel.resize(2, nSides); 

    for (int ie = 0; ie < nElements; ie++) {
        
        for (int j = 0; j < i34[ie]; j++) {
            const int jsj = elside(j, ie);  

            // const int n1 = elnode(nxQ(0, j, i34[ie]), ie);
            // const int n2 = elnode(nxQ(1, j, i34[ie]), ie);

            if (ic3(j,ie)==0 || (ic3(j,ie)>0)) { 
            // // && iegb<ielg(max(1,ic3(j,ie)))) ||
            //     // (ic3(j,ie)<0 && iegb<iabs(ic3(j,ie)))) {

                isdel(0,jsj)=ie;
                isdel(1,jsj)=ic3(j,ie);
            //     // isidenode(1,jsj)=n1
            //     // isidenode(2,jsj)=n2
            } else {
                // old augmented subdomain boundary side (reverse orientation)

                isdel(0,jsj)=ic3(j,ie);
                isdel(1,jsj)=ie;
            //     // isidenode(1,jsj)=n2
            //     // isidenode(2,jsj)=n1
            }
        }
    }


    // for (int is = 0; is < nSides; is++) {
    //     isdel(0, is) = js(0, is);
    //     isdel(1, is) = ic3(0, is);
    // }

    

    print("");

}


// void buildSideData(int nElements, int nSides, 
//                   const Array2D<int>& elnode, const Array2D<int>& elside, const std::vector<int>& i34,
//                   const Array2D<int>& ic3, const ArrayOfRows<int>& vertexElements,
//                   const Array2D<int>& isdel, const Array2D<int>& isidenode,
//                   const Array2D<double>& xcj, const Array2D<double>& ycj, const Array2D<double>& zcj,
//                   const Array2D<double>& dps, const Array2D<double>& side) {
    
//     for (int ie = 0; ie < nElements; ie++) {
        
//         for (int j = 0; j < i34[ie]; j++) {
//             int jsj = elside(j, ie);  // Side index
//             int n1 = elnode(nxQ(0, j, i34[ie]), ie);  // First node of side
//             int n2 = elnode(nxQ(1, j, i34[ie]), ie);  // Second node of side
            
//             // Check if this is a new local side
//             if (ic3(j, ie) == 0 || 
//                (ic3(j, ie) > 0 && iegb < ielg(std::max(0, ic3(j, ie)))) ||
//                (ic3(j, ie) < 0 && iegb < std::abs(ic3(j, ie)))) {
                
//                 // Store element indices and nodes in original orientation
//                 isdel[0, jsj] = ie;
//                 isdel[1, jsj] = ic3(j, ie);
//                 isidenode[0, jsj] = n1;
//                 isidenode[1, jsj] = n2;
//             } else {
//                 // Old augmented subdomain boundary side (reverse orientation)
//                 isdel[0, jsj] = ic3(j, ie);
//                 isdel[1, jsj] = ie;
//                 isidenode[0, jsj] = n2;  // Note reversed order
//                 isidenode[1, jsj] = n1;
//             }
            
//             // Calculate side center coordinates
//             xcj[jsj] = (xnd[n1] + xnd[n2]) * 0.5;
//             ycj[jsj] = (ynd[n1] + ynd[n2]) * 0.5;
//             if (ics == 2) {  // 3D case 
//                 zcj[jsj] = (znd[n1] + znd[n2]) * 0.5;
//             }
            
//             // Calculate depth and length
//             dps[jsj] = (dp[n1] + dp[n2]) * 0.5;
            
//             // Calculate side length
//             double dx = xnd[n2] - xnd[n1];
//             double dy = ynd[n2] - ynd[n1];
//             double dz = (ics == 2) ? znd[n2] - znd[n1] : 0.0;
//             side[jsj].length = std::sqrt(dx*dx + dy*dy + dz*dz);
            
//             if (side[jsj].length == 0.0) {
//                 throw std::runtime_error("Zero length side: " + std::to_string(jsj));
//             }
//         }
//     }
// }


void reorderElementsCounterClockwise(int np, const std::vector<int>& isbnd, 
                                   ArrayOfRows<int>& vertexElements,
                                   const Array2D<int>& elnode, const std::vector<int>& i34,
                                   const Array2D<int>& ic3) {
    for (int i = 0; i < np; i++) {
        progressBar("Reordering node elements counter-clockwise", i, np, ansi::LIGHT_BLUE, 20);
        // If it is a boundary node, starts with a boundary element
        if (isbnd[i] != -1) {
            // Find first boundary element and place it first
            int icount = 0;
            for (size_t ielocal = 0; ielocal < vertexElements.rowSize(i); ielocal++) {
                const int ie = vertexElements(i, ielocal);

                // Find local index of current node i in element ie
                int ilocal = -1;
                for (int l = 0; l < i34[ie]; l++) {
                    if (elnode(l, ie) == i) {
                        ilocal = l;
                        break;
                    }
                }
                if (ilocal == -1) {
                    throw std::runtime_error("Node " + std::to_string(i+1) + 
                                          " not found in its own element " + 
                                          std::to_string(ie+1));
                }

                // If ilocal is a boundary node and ilocal-1 is too, then this is a boundary element
                if (ic3(nxQ(-1, ilocal, i34[ie]), ie) == -1) {
                    icount++;
                    vertexElements(i, 0) = ie;
                }
            }

            if (icount != 1) {
                throw std::runtime_error("Illegal boundary node " + std::to_string(i+1) + 
                                       ": found " + std::to_string(icount) + 
                                       " boundary elements");
            }
        }

        // For all nodes: recompute remaining elements in counter-clockwise order
        for (size_t j = 1; j < vertexElements.rowSize(i); j++) {
            // Get previous element
            const int ie = vertexElements(i, j - 1);

            // Find local index of node i in element ie
            int ilocal = -1;
            for (int l = 0; l < i34[ie]; l++) {
                if (elnode(l, ie) == i) {
                    ilocal = l;
                    break;
                }
            }
            if (ilocal == -1) {
                throw std::runtime_error("Node " + std::to_string(i+1) + 
                                         " not found in its own element " + 
                                         std::to_string(ie+1));
            }

            // Find next element by going two sides back counter-clockwise
            const int newElement = ic3(nxQ(i34[ie] - 3, ilocal, i34[ie]), ie);
            if (newElement == -1) {
                throw std::runtime_error("Incomplete ball around node " + 
                                         std::to_string(i+1) + 
                                         " at surrounding local element " + 
                                         std::to_string(j));
            }

            vertexElements(i, j) = newElement;
        }
    }
    print("");
}

void reorderPartitionNodes(
    const std::vector<idx_t>& partitions,
    const AdjacencyList& fullAdjList,
    int targetPartition) 
{

        // Set up METIS options specifically for node reordering
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    
    // Optimize for fill reduction
    options[METIS_OPTION_NSEPS] = 1;      // Number of separators: fewer separators = faster
    options[METIS_OPTION_NITER] = 5;      // Number of refinement iterations
    options[METIS_OPTION_UFACTOR] = 1;    // Ultra-tight balance for better locality
    options[METIS_OPTION_COMPRESS] = 1;   // Compress graph for faster processing
    options[METIS_OPTION_CCORDER] = 1;    // Connected components ordered separately

    // Count nodes in this partition
    int nNodesInPartition = std::count(partitions.begin(), partitions.end(), targetPartition);
    
    // Create mapping of old to new indices
    std::vector<int> localToGlobal;
    localToGlobal.reserve(nNodesInPartition);
    
    // Create reverse mapping (global to local)
    std::vector<int> globalToLocal(partitions.size(), -1);
    int localIdx = 0;
    for (size_t i = 0; i < partitions.size(); i++) {
        if (partitions[i] == targetPartition) {
            localToGlobal.push_back(i);
            globalToLocal[i] = localIdx++;
        }
    }

    // Build local adjacency list
    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    
    xadj.push_back(0);
    for (const auto& globalIdx : localToGlobal) {
        // For each node in partition
        for (int j = fullAdjList.xadj[globalIdx]; j < fullAdjList.xadj[globalIdx + 1]; j++) {
            int neighbor = fullAdjList.adjncy[j];
            // Only include edges to nodes in same partition
            if (partitions[neighbor] == targetPartition) {
                adjncy.push_back(globalToLocal[neighbor]);
            }
        }
        xadj.push_back(adjncy.size());
    }

    // Prepare arrays for METIS
    std::vector<idx_t> perm(nNodesInPartition);    // Will hold new ordering
    std::vector<idx_t> iperm(nNodesInPartition);   // Inverse permutation

#ifdef popo
    // Call METIS
    printFormatted(ansi::yellow("Reordering partition %d with %d nodes"), targetPartition, nNodesInPartition);
    int ret = METIS_NodeND(
        &nNodesInPartition,  // number of vertices
        xadj.data(),         // adjacency structure
        adjncy.data(),       // adjacency structure
        nullptr,             // weights of vertices
        options,             // options
        &perm[0],            // permutation array
        &iperm[0]            // inverse permutation array
    );

    if (ret != METIS_OK) {
        throw std::runtime_error("METIS_NodeND failed");
    }

    for (int i = 0; i < 100; i++) {
        printFormatted(ansi::green("perm[%d] = %d, iperm[%d] = %d"), i, perm[i], i, iperm[i]);
    }
#endif
    // perm now contains the new ordering
    // perm[old_local_idx] = new_local_idx
    // iperm[new_local_idx] = old_local_idx
}

void partitionNodesAndSides(int nElements, int nNodes, int nSides,
                            const Array2D<int>& elnode, const std::vector<int>& i34, 
                            const Array2D<int>& ic3, const Array2D<int>& isdel,
                            const std::vector<idx_t>& partitions,
                            const ArrayOfRows<int>& vertexElements,
                            std::vector<int> &partNodes, std::vector<int> &partSides)
{
    // Assigns nodes and sides also to different partitions. This is done by choosing
    // the partition of one of the elements that the node belongs to.

    partNodes.resize(nNodes);
    partSides.resize(nSides);


    // Assign nodes to partitions
    for (int in = 0; in < nNodes; in++) {
        const int numElementsPerNode = vertexElements.rowSize(in);

        // Choses one of the elements that the node belongs to and assigns the partition
        // of that element to the node.
        const int le = rand() % numElementsPerNode;

        const int ie = vertexElements(in, le);
        partNodes[in] = partitions[ie];
    }

    printFormatted(ansi::blue("Partitfdfioning sides = %d"), nSides);

    // Assign sides to partitions
    for (int is = 0; is < nSides; is++) {

        const int le = rand() % 2;
        partSides[is] = partitions[isdel(le, is)];
    }

}


// Helper functions to convert C++ arrays to Python arrays
template<typename T>
py::array_t<T> toPyArray(const std::vector<T>& data) {
    auto array = py::array_t<T>(data.size());
    auto ptr = static_cast<T*>(array.request().ptr);
    std::copy(data.begin(), data.end(), ptr);
    return array;
}

template<typename T>
py::array_t<T> toPyArray(const Array2D<T>& data) {
    auto array = py::array_t<T>({data.rows(), data.cols()});
    auto ptr = static_cast<T*>(array.request().ptr);
    std::copy(&data(0,0), &data(0,0) + data.rows() * data.cols(), ptr);
    return array;
}

std::tuple<py::array_t<int32_t>,    // i34
           py::array_t<int32_t>,    // partElements
           py::array_t<int32_t>,    // partNodes
           py::array_t<int32_t>,    // partSides
           py::array_t<int32_t>,    // elnode
           py::array_t<int32_t>,    // elside
           py::array_t<double>,     // nodes
           py::array_t<int32_t>,    // ic3
           py::array_t<int32_t>>    // isdel
readMesh(const std::string& filename, idx_t nPartitions) {


    // Read the header information
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Could not open file");
    
    // Read mesh name
    std::string meshName;
    std::getline(file, meshName);

    // Read element and node counts
    int nElements, nNodes, nSides;

    auto numElementsNodes = readInts(file);
    nElements = numElementsNodes[0];
    nNodes    = numElementsNodes[1];


    // Initialize arrays for node and element data
    Array2D<double> nodes(nNodes, 3);      // [x, y, depth]
    std::vector<int> i34(nElements);    // Number of nodes per element
    Array2D<int> elnode(4, nElements);  // Node indices for each element
    Array2D<int> elside(4, nElements);

    printFormatted("%s: %i, %s: %i", ansi::red("n_nodes: ").c_str(), nNodes, ansi::green("n_elements: ").c_str(), nElements); 

    // Read node information
    std::string line;
    for (int i = 0; i < nNodes; i++) {
        progressBar("Reading nodes", i, nNodes, ansi::RED, 20);
        std::getline(file, line);
        std::istringstream nodeIss(line); 
        int index;
        nodeIss >> index;  

        try {
            for (int j = 0; j < 3; j++) {
                nodeIss >> nodes(i, j);
            }
        } catch (const std::exception& e) {
            std::cerr << "Error reading node " << i << std::endl;
            std::cerr << "Line: " << line << std::endl;
            throw;
        }
    }

    // Read element information
    print("");
    for (int i = 0; i < nElements; i++) {
        progressBar("Reading elements", i, nElements, ansi::GREEN, 20);
        std::getline(file, line);
        std::istringstream elemIss(line);
        int index;  // Skip the index
        elemIss >> index;
        try {
            elemIss >> i34[i];  // Number of nodes in element
            for (int j = 0; j < i34[i]; j++) {
                int nodeIndex;
                elemIss >> nodeIndex;
                elnode(j, i) = nodeIndex - 1;  // Convert to 0-based indexing
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error reading element " << i << std::endl;
            std::cerr << "Line: " << line << std::endl;
            throw;
        }
    }
    print("");

    // Read open boundary information
    // Read number of open boundary segments and total boundary nodes
    int nope = readInt(file);  
    int neta = readInt(file);  

    // First pass: count max nodes per segment and verify total
    int mnond = 0;  // max number of nodes per segment
    int nt = 0;     // total node count
    std::vector<std::vector<int>> segmentNodes;  // temporarily store nodes per segment

    // Process open boundary segments
    for (int k = 0; k < nope; k++) {
        progressBar("Reading open boundaries", k, nope, "\033[33m", 20);

        int nn = readInt(file);  // number of nodes in this segment
        mnond = std::max(mnond, nn);
        nt += nn;
        // Read node indices for this segment
        std::vector<int> segment(nn);
        for (int i = 0; i < nn; i++) {
            segment[i] = readInt(file) - 1;  // Convert to 0-based indexing
        }
        segmentNodes.push_back(segment);
        

    }
    print("");
    
    if (neta != nt) {
        char msg[256];
        snprintf(msg, sizeof(msg), "neta (%d) != total # of open boundary nodes (%d)", neta, nt);
        throw std::runtime_error(msg);
    }

    // Allocate arrays for open boundary segments
    std::vector<int> nond(nope, 0);                    // nodes per segment
    Array2D<int> iond(nope, mnond);                    // node indices for each segment
    std::vector<int> isbnd(nNodes, -1);                // which segment each node belongs to

    // Fill the arrays using stored segment_nodes
    for (size_t k = 0; k < segmentNodes.size(); k++) {
        const auto& nodes = segmentNodes[k];
        int nn = nodes.size();
        nond[k] = nn;
        
        // Copy node indices
        for (int i = 0; i < nn; i++) {
            iond(k, i) = nodes[i];
        }

        // Mark boundary nodes
        for (int ip : nodes) {
            isbnd[ip] = k;
        }

        // Check for looped boundaries
        if (nodes[0] == nodes[nodes.size()-1]) {
            char msg[256];
            snprintf(msg, sizeof(msg), "Looped open boundary: %lud", k+1);
            throw std::runtime_error(msg);
        }
    }

    if (nope > 0) {
        std::string msg = "Read " + std::to_string(nope) + " open boundary segments with " + std::to_string(nt) + " total nodes";
        print(msg);
    }

    // Read land boundary information
    // Read number of land boundary segments and total boundary nodes
    std::getline(file, line);
    std::istringstream iss1(line);
    int nland;
    iss1 >> nland;

    std::getline(file, line);
    std::istringstream iss2(line);
    int nvel;
    iss2 >> nvel;

    // First pass: count max nodes per segment and verify total
    int mnlnd = 0;  // max number of nodes per segment
    nt = 0;     // total node count reset for land boundaries
    segmentNodes.clear();  // temporarily store nodes per segment

    for (int k = 0; k < nland; k++) {
        progressBar("Reading land boundaries", k, nland, ansi::BLUE, 20);
        // Read number of nodes in this segment
        std::getline(file, line);
        std::istringstream iss(line);
        int nn;
        iss >> nn;
        mnlnd = std::max(mnlnd, nn);
        nt += nn;

        // Read and store node indices
        std::vector<int> nodes(nn);
        for (int i = 0; i < nn; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> nodes[i];
            nodes[i]--;  // Convert to 0-based indexing
        }
        segmentNodes.push_back(nodes);

 
    }

    print("");
    if (nvel != nt) {
        throw std::runtime_error("nvel (" + std::to_string(nvel) + 
                               ") != total # of land boundary nodes (" + 
                               std::to_string(nt) + ")");
    }

    // Allocate arrays for land boundary segments
    std::vector<int> nlnd(nland, 0);      // nodes per segment
    std::vector<std::vector<int>> ilnd(nland, std::vector<int>(mnlnd, 0));  // node indices for each segment

    // Fill the arrays using stored segment_nodes
    for (size_t k = 0; k < segmentNodes.size(); k++) {
        const auto& nodes = segmentNodes[k];
        int nn = nodes.size();
        nlnd[k] = nn;
        
        // Copy node indices
        for (int i = 0; i < nn; i++) {
            ilnd[k][i] = nodes[i];
        }

        // Mark boundary nodes and check for overlap with open boundaries
        for (int ip : nodes) {
            if (isbnd[ip] == -1) {
                isbnd[ip] = -2;  // overlap with open boundary
            }
        }
    }

    printFormatted(ansi::blue("Read %d land boundary segments with %d total nodes"), nland, nt);
    print("Done reading mesh.");


    // Process the mesh
    ArrayOfRows<int> vertexElements;


    buildVertexElements(elnode, i34, nElements, nNodes, vertexElements);

    AdjacencyList adjList = buildAdjacencyList(elnode, i34, vertexElements, nElements);

    Array2D<int> js(4, 1, -1);
    Array2D<int> ic3(4, nElements, -1);
    Array2D<int> isdel(2, 1, -1); // We still don't know the number of sides

    buildElementSideElementTable(elnode, i34, vertexElements, nElements, nNodes, 
                                 elside, isdel, ic3, nSides);


    reorderElementsCounterClockwise(nNodes, isbnd, vertexElements, elnode, i34, ic3);


    

    printFormatted(ansi::green("Partitioning mesh using METIS into %d partitions, nElements = %d"), nPartitions, nElements);

    // Prepare METIS arrays
    idx_t nvtxs = nElements;  // number of vertices (elements in our case)
    idx_t ncon = 1;          // number of balancing constraints

    std::vector<idx_t> partElements(nvtxs);  // array to store partition numbers
    
    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    
    printFormatted(ansi::red("Partitioning mesh using METIS into %d partitions"), nPartitions);
    // Call METIS directly with adjList's data
    idx_t objval;
    int ret = METIS_PartGraphKway(
        &nvtxs,
        &ncon,
        adjList.xadj.data(),     // Use directly from adjList
        adjList.adjncy.data(),   // Use directly from adjList
        nullptr,
        nullptr,
        adjList.adjwgt.data(),   // Use directly from adjList
        &nPartitions,
        nullptr,
        nullptr,
        options,
        &objval,
        partElements.data()
    );

    if (ret != METIS_OK) {
        throw std::runtime_error("METIS partitioning failed");
    }

    printFormatted(ansi::green("METIS partitioning completed with Edge cut ratio: %f"), objval / float(adjList.ntedge));

    bool reorderPartitions = false;
    if (reorderPartitions) {
        printFormatted(ansi::yellow("Reordering %d partitions"), nPartitions);
        for (int ipart = 0; ipart < nPartitions; ipart++) {
            progressBar("Reordering partition nodes", ipart, nPartitions, ansi::YELLOW, 20);
            reorderPartitionNodes(partElements, adjList, ipart);
        }
    }

    // Assigns nodes and sides also to different partitions.
    std::vector<int> partNodes;
    std::vector<int> partSides;

    partitionNodesAndSides(nElements, nNodes, nSides, elnode, i34, ic3, isdel, partElements, vertexElements, partNodes, partSides);

    
    return std::make_tuple(
        toPyArray(i34),
        toPyArray(partElements),
        toPyArray(partNodes),
        toPyArray(partSides),
        toPyArray(elnode),
        toPyArray(elside),
        toPyArray(nodes),
        toPyArray(ic3),
        toPyArray(isdel)

    );
}

PYBIND11_MODULE(test, m) {
    m.def("readMesh", &readMesh, "Read a mesh file");
}
