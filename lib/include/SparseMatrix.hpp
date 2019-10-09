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
    T* link;    
};

template<typename T>
class SparseMatrix {
    unsigned height, width;
    Node<T> ** data;
public:
    // default constructor
    SparseMatrix():height(0), width(0), data(nullptr) {}
    //Copy constructor
    SparseMatrix(const SparseMatrix &sm);
    //Move constructor
    SparseMatrix(SparseMatrix &&sm);

    SparseMatrix(const unsigned &height, const unsigned &width, const vector<vector<Node<T>>> &data);
    

     // Initilize from dense matrix located in file
    static SparseMatrix<T> fromFileDense(string file_name);
    static SparseMatrix<T> fromFileCSC(string file_name);

    const vector<Node<T>>& getVector(const unsigned &index) const { return data[index]; }
    T& operator()(const unsigned &i, const unsigned &j) { return data[i][j].value; }
    const T& operator()(const unsigned &i, const unsigned &j) const { return data[i][j].value; }
    unsigned getIndex(const unsigned &i, const unsigned &j) { return data[i][j].index; }
    unsigned getIndex(const unsigned &i, const unsigned &j) const { return data[i][j].index; }
    T* getAddress(const unsigned &i, const unsigned &j) { return &data[i][j].value; }
    T* getLink(const unsigned &i, const unsigned &j) { return data[i][j].link; }
    T* getLink(const unsigned &i, const unsigned &j) const { return data[i][j].link; }
    const unsigned getHeight() const { return height; }
    const unsigned getWidth() const { return width; }
    const unsigned getRowSize(const unsigned &i) const { return data[i].size(); }
    
    void clear();
    void printSparse(const unsigned &spacing = 3) const;
    void printDense(const unsigned &spacing = 4) const;
    void insertRow(const vector<Node<T>> row, const unsigned &index) {
        data.insert(data.begin() + index, row);
        height++;
    }
    void swapRows(const unsigned& i, const unsigned& j){
        data[i].swap(data[j]);
    }
    void link(const unsigned &i, const unsigned &j, T* address) { data[i][j].link = address; }
    
    SparseMatrix getTransposed() const;
    SparseMatrix getSquared() const;
    SparseMatrix getInverse() const;

};

#endif // SPARSE_MATRIX_H