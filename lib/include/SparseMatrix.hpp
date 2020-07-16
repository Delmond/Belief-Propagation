#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <string>
#include <algorithm>

using namespace std;

template<typename T>
struct Node {
    T value;
    unsigned index;
};

template<typename T>
class SparseMatrix {
    unsigned height, width, nnz;
    T *values;
    unsigned *indices, *offset, *ptr;
public:
    // default constructor
    SparseMatrix():height(0), width(0), nnz(0), values(nullptr), indices(nullptr), offset(nullptr), ptr(nullptr) {}
    //Copy constructor
    SparseMatrix(const SparseMatrix &sm);
    //Move constructor
    SparseMatrix(SparseMatrix &&sm);

    SparseMatrix(const unsigned &height, const unsigned &width, const vector<vector<Node<T>>> &data);
    
    SparseMatrix<T>& operator=(const SparseMatrix<T>&) = delete;
    SparseMatrix<T>& operator=(SparseMatrix<T>&&);
     // Initilize from dense matrix located in file
    static SparseMatrix<T> fromFileDense(string file_name);
    static SparseMatrix<T> fromFileCSC(string file_name);
    T& operator()(const unsigned &i, const unsigned &j) { return values[offset[i] + j]; }
    const T& operator()(const unsigned &i, const unsigned &j) const { return values[offset[i] + j]; }
    unsigned getDirectIndex(const unsigned &i, const unsigned &j) { return offset[i] + j; }
    unsigned getDirectIndex(const unsigned &i, const unsigned &j) const { return offset[i] + j; }
    
    unsigned getIndex(const unsigned &i, const unsigned &j) { return indices[offset[i] + j]; }
    unsigned getIndex(const unsigned &i, const unsigned &j) const { return indices[offset[i] + j]; }
    
    T& operator[](const unsigned &i){ return values[i]; }
    const T& operator[](const unsigned &i) const { return values[i]; }


    const unsigned getHeight() const { return height; }
    const unsigned getWidth() const { return width; }
    const unsigned getRowSize(const unsigned &i) const { return offset[i + 1] - offset[i]; }
    const unsigned getNNZ() const { return nnz; }     
    unsigned& getLink(const unsigned &i, const unsigned &j) { return ptr[offset[i] + j]; }
    const unsigned& getLink(const unsigned &i, const unsigned &j) const { return ptr[offset[i] + j]; }
    unsigned& getDirectLink(const unsigned &i) { return ptr[i]; }
    const unsigned& getDirectLink(const unsigned &i) const { return ptr[i]; }

    void clear();
    void printSparse(const unsigned &spacing = 3) const;
    void printDense(const unsigned &spacing = 4, const unsigned &precision = 3) const;

    T* getVector(const unsigned &i) { return values + offset[i]; }
    // void insertRow(const vector<Node<T>> row, const unsigned &index) {
    //     data.insert(data.begin() + index, row);
    //     height++;
    // }
    // void link(const unsigned &i, const unsigned &j, T* address) { data[i][j].link = address; }
    vector<vector<Node<T>>> getData() const;
    SparseMatrix getTransposed() const;
    SparseMatrix getSquared() const;
    SparseMatrix getInverse() const;

    ~SparseMatrix();
};

#endif // SPARSE_MATRIX_H