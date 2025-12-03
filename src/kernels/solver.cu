#include <cstdint>
#include "meshdata.h"

/**
 * Sparse matrix in CRS format
 * 
 * Conceptually, this class represents the global sparse matrix, but only the local rows are stored in local memory.
 * The global matrix is stored in a compressed row storage (CRS) format, where the non-zero values are stored in a single array,
 * and the column indices and row starts are stored in separate arrays.
 */
struct SparseMatrix {

    const Mesh &mesh;

    
    // Local memory only
    double*  values;           // Non-zero values for local rows
    int32_t* columnIndices;   // Global column indices for local values
    int32_t* rowStarts;       // Start indices for local rows
    
    // Counts
    const int32_t localRows;       // Number of rows in this GPU (same as local nodes)
    const int32_t localNonZeros;   // Number of non-zeros in this GPU

    SparseMatrix(const Mesh &mesh) : mesh(mesh) {
        // Initialize the sparse matrix structure

        // Count the number of local rows and non-zeros
        localRows = mesh.np;
        localNonZeros = 0;
        for (int32_t i = 0; i < mesh.np; i++) {
            const Node &node = mesh.nodes[i];
            localNonZeros += node.numNodes();
        }

        // Initialize the local memory arrays
        values        = (double* )malloc(localNonZeros   * sizeof(double));
        columnIndices = (int32_t*)malloc(localNonZeros   * sizeof(int32_t));
        rowStarts     = (int32_t*)malloc((localRows + 1) * sizeof(int32_t));
    }

    ~SparseMatrix() {
        free(values);
        free(columnIndices);
        free(rowStarts);
    }

    void mult(const double* x, double* y) const
    {
        // TODO: Implement the matrix-vector multiplication
    }
};