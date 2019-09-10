#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H


#include <vector>
#include <string>

using namespace std;

template<typename T>
struct Node {
    T value;
    unsigned index;
    T* link;    
};

template<typename T>
class SparseMatrix {
    unsigned height, width;
    vector<vector<Node<T>>> data;
public:
    // default constructor
    SparseMatrix():height(0), width(0), data() {}
    // Initilize from dense matrix located in file
    SparseMatrix(string file_name);
    //Copy constructor
    SparseMatrix(SparseMatrix &sm);
    //Move constructor
    SparseMatrix(SparseMatrix &&sm);

    const vector<Node<T>>& getVector(const unsigned &index) const;
    T& operator()(const unsigned &i, const unsigned &j) { return data[i][j].value; }
    const unsigned getHeight() const { return height; }
    const unsigned getWidth() const { return width; }
    const unsigned getRowSize(const unsigned &i) const { return data[i].size(); }
    void printSparse(const unsigned &spacing = 3) const;
    void printDense(const unsigned &spacing = 4) const;
    unsigned getIndex(const unsigned &i, const unsigned &j) { return data[i][j].index; }
    T* getAddress(const unsigned &i, const unsigned &j) { return &data[i][j].value; }
    void link(const unsigned &i, const unsigned &j, T* address) { data[i][j].link = address; }
    
    SparseMatrix getTransposed() const;

};

#endif // SPARSE_MATRIX_H