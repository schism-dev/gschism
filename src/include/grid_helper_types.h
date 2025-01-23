// Data structures to store grid related data.
// Since C++ can't have 2D arrays we use parenthetical indexing.

#ifdef DEBUG // Define DEBUG in your debug build configuration
#define ARRAY_BOUNDS_CHECK(i, j, sizeI, sizeJ) \
    if (i >= sizeI || j >= sizeJ || i < 0 || j < 0) \
        throw std::out_of_range("Array2D index out of bounds: (" + std::to_string(i) + ", " + std::to_string(j) + ")");
#else
#define ARRAY_BOUNDS_CHECK(i, j, sizeI, sizeJ) // No checks in release mode
#endif



struct AdjacencyList {
    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    std::vector<idx_t> adjwgt;
    idx_t ntedge;

    ~AdjacencyList() {
        xadj.clear();
        adjncy.clear();
        adjwgt.clear();
    }
};


// 2D "rectangular array class.
template <typename T>
class Array2D {
private:
    T* data;
    size_t sizeI;
    size_t sizeJ;

public:
    Array2D(size_t sizeI, size_t sizeJ, T fill = 0) : sizeI(sizeI), sizeJ(sizeJ) {
        data = (T*)malloc(sizeI * sizeJ * sizeof(T));
        std::fill(data, data + sizeI * sizeJ, fill);
    } 

    void resize(size_t sizeI, size_t sizeJ, T fill = 0) {
        free(data);

        data = (T*)malloc(sizeI * sizeJ * sizeof(T));
        std::fill(data, data + sizeI * sizeJ, fill);

        this->sizeI = sizeI;
        this->sizeJ = sizeJ; 
    }

    ~Array2D() {
        free(data);
    } 

    // Access operator
    T& operator()(size_t i, size_t j) {
        ARRAY_BOUNDS_CHECK(i, j, sizeI, sizeJ);

        return data[i * sizeJ + j];
    }

    // Const access operator
    const T& operator()(size_t i, size_t j) const {
       ARRAY_BOUNDS_CHECK(i, j, sizeI, sizeJ);

        return data[i * sizeJ + j];
    }

    // Get dimensions
    size_t rows() const { return sizeI; }
    size_t cols() const { return sizeJ; }
    
    // Get raw data pointer
    T* data_ptr() { return data; }
    const T* data_ptr() const { return data; }
};

// 2D array with arbitrary length rows class.
template <typename T>
class ArrayOfRows {
private:
    T* data;
    size_t sizeI;
    size_t* rowStarts;

public:
    ArrayOfRows(size_t sizeI, size_t* rowLengths) : sizeI(sizeI) {
        rowStarts = (size_t*)malloc((sizeI + 1) * sizeof(size_t));

        // Compute the start index for each row
        rowStarts[0] = 0;
        for (size_t i = 1; i <= sizeI; i++) {
            rowStarts[i] = rowStarts[i-1] + rowLengths[i-1];
        }
        data = (T*)malloc(rowStarts[sizeI] * sizeof(T));
    }

    ArrayOfRows(const std::vector<size_t>& rowLengths) : ArrayOfRows(rowLengths.size(), rowLengths.data()) {}

    ArrayOfRows(void): data(nullptr), sizeI(0), rowStarts(nullptr) {}

    void rebuild(const std::vector<size_t>& rowLengths) {
      // Rebuild the data structure for a new set of row lengths. It clears the old data.
        free(data);
        free(rowStarts);
        sizeI = rowLengths.size();
        rowStarts = (size_t*)malloc((sizeI + 1) * sizeof(size_t));
        rowStarts[0] = 0;
        for (size_t i = 1; i <= sizeI; i++) {
            rowStarts[i] = rowStarts[i-1] + rowLengths[i-1];
        }
        data = (T*)malloc(rowStarts[sizeI] * sizeof(T));
    }

    ~ArrayOfRows() {
        free(rowStarts);
        free(data);
    }

    // Access operator
    T& operator()(size_t i, size_t j) {
        assert(i < sizeI && j < rowLengths[i] && i >= 0 && j >= 0);
        return data[rowStarts[i] + j];
    }

    // Const access operator
    const T& operator()(size_t i, size_t j) const {
        assert(i < sizeI && j < rowLengths[i] && i >= 0 && j >= 0);
        return data[rowStarts[i] + j];
    }

    size_t rowSize(size_t i) const { return rowStarts[i+1] - rowStarts[i]; }
};
