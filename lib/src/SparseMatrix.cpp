#include "../include/SparseMatrix.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>

constexpr double EPSILON = 1e-7;

using namespace std;

template<typename T>
SparseMatrix<T> SparseMatrix<T>::fromFileDense(string file_name) {
    unsigned height, width;
    vector<vector<Node<T>>> data;

    std::ifstream fileInput(file_name);
    if(fileInput.fail()){
        throw domain_error("Failed opening file in Sparse Matrix fromFileDense static method.");
    }

    fileInput >> height >> width;
    data.resize(height);
    T value;


    for(unsigned i = 0; i < height; i++){    
        for(unsigned j = 0; j < width; j++){
            fileInput >> value;
            if(abs(value) < EPSILON)
                continue;

            data[i].push_back({value, j, nullptr});
        }
    }
    return SparseMatrix<T>(height, width, data);
}

template<typename T>
SparseMatrix<T> SparseMatrix<T>::fromFileCSC(string file_name){
    unsigned height(0), width(0);
    vector<vector<Node<T>>> data;

    ifstream fileInput(file_name);
    if(fileInput.fail()){
        throw domain_error("Failed opening file in Sparse Matrix fromFileCSC static method.");
    }
    unsigned row, column;
    T value;
    while(true){
        fileInput >> row >> column >> value;
        if(fileInput.eof())
            break;
        // for 1-Indexed
        column--;
        row--;
        
        if(width <= column) width = column + 1;
        if(data.size() <= row) data.resize(row + 1);

        data[row].push_back({value, column, nullptr});
    }
    height = data.size();
    return SparseMatrix<T>(height, width, data);
}

template<typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix &sm):
    height(sm.height),
    width(sm.width),
    data(sm.data)
    {}


template<typename T>
SparseMatrix<T>::SparseMatrix(SparseMatrix<T> &&sm):
    height(sm.height),
    width(sm.width),
    data(move(sm.data))
    {}




template<typename T>
void SparseMatrix<T>::clear(){
    for(unsigned i = 0; i < data.size(); i++){
        for(unsigned j = 0; j < data[i].size(); j++){
            data[i][j].value = T();
        }
    }

}

template<typename T>
void SparseMatrix<T>::printSparse(const unsigned &spacing) const {
    for(unsigned i = 0; i < data.size(); i++){
        for (unsigned j = 0; j < data[i].size(); j++){
            cout << setw(spacing) <<data[i][j].value << " (" << i << ", " << data[i][j].index << ")" << endl; 
        }
    }
}

template<typename T>
void SparseMatrix<T>::printDense(const unsigned &spacing) const {
    unsigned position;
    for(unsigned i = 0; i < data.size(); i++){
        position = 0;
        for(unsigned j = 0; j < data[i].size(); j++){
            while(position != data[i][j].index){
                cout << setw(spacing) << 0;
                position++;
            }
            cout << setw(spacing) << data[i][j].value;
            position++;
        }
        while(position != width){
            cout << setw(spacing) << 0;
            position++;
        }
        cout << "\n";
    }
}

template<typename T>
SparseMatrix<T> SparseMatrix<T>::getTransposed() const {
    SparseMatrix<T> transposed;
    transposed.height = width;
    transposed.width = height;
    transposed.data.resize(width);

    for(unsigned i = 0; i < height; i++){
        for(unsigned j = 0; j < data[i].size(); j++){
            unsigned column = data[i][j].index;
            transposed.data[column].push_back({data[i][j].value, i, nullptr});
        }
    }
    return transposed;
}
template<typename T>
SparseMatrix<T> SparseMatrix<T>::getSquared() const {
    SparseMatrix<T> squared(*this);
    for(unsigned i = 0; i < squared.getHeight(); i++){
        for(unsigned j = 0; j < squared.getRowSize(i); j++){
            squared(i, j) *= squared(i, j);
        }
    }
    return squared;
}
template<typename T>
SparseMatrix<T> SparseMatrix<T>::getInverse() const {
    SparseMatrix<T> inverse(*this);
    for(unsigned i = 0; i < inverse.getHeight(); i++){
        for(unsigned j = 0; j < inverse.getRowSize(i); j++){
            inverse(i, j) = 1.0/inverse(i, j);
        }
    }
    return inverse;
}

template class SparseMatrix<float>;
template class SparseMatrix<double>;