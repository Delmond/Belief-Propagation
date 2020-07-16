#include "../include/SparseMatrix.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <unistd.h>

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

            data[i].push_back({value, j});
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

        data[row].push_back({value, column});
    }
    height = data.size();
    return SparseMatrix<T>(height, width, data);
}

template<typename T>
SparseMatrix<T>::SparseMatrix(const unsigned &height, const unsigned &width, const vector<vector<Node<T>>> &data)
:height(height), width(width){
    nnz = 0;
    for(unsigned int i = 0; i < data.size(); i++){
        nnz += data[i].size();
    }
    values = new T[nnz];
    indices = new unsigned[nnz];
    ptr = new unsigned[nnz];
    offset = new unsigned[height + 1];
    unsigned index = 0;
    offset[0] = 0;
    for(unsigned int i = 0; i < data.size(); i++){
        for(unsigned int j = 0; j < data[i].size(); j++){
            values[index] = data[i][j].value;
            indices[index] = data[i][j].index;
            index++;
        }
        offset[i + 1] = offset[i] + data[i].size();
    }

}

template<typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix &sm):
    height(sm.height),
    width(sm.width),
    nnz(sm.nnz)
    {
        values = new T[nnz];
        indices = new unsigned[nnz];
        ptr = new unsigned[nnz];
        offset = new unsigned[height + 1];
        copy(sm.values, sm.values + nnz, values);
        copy(sm.indices, sm.indices + nnz, indices);
        copy(sm.ptr, sm.ptr + nnz, ptr);
        copy(sm.offset, sm.offset + height + 1, offset);
    }


template<typename T>
SparseMatrix<T>::SparseMatrix(SparseMatrix<T> &&sm):
    height(sm.height),
    width(sm.width),
    nnz(sm.nnz),
    values(sm.values),
    indices(sm.indices),
    offset(sm.offset),
    ptr(sm.ptr)
    {
        sm.values = nullptr;
        sm.indices = nullptr;
        sm.offset = nullptr;
        sm.ptr = nullptr;
    }
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(SparseMatrix<T>&& sm){

    height = sm.height;
    width = sm.width;
    nnz = sm.nnz;
    values = sm.values;
    indices = sm.indices;
    offset = sm.offset;
    ptr = sm.ptr;

    sm.values = nullptr;
    sm.indices = nullptr;
    sm.offset = nullptr;
    sm.ptr = nullptr;

    return *this;

}

template<typename T>
SparseMatrix<T>::~SparseMatrix(){
    delete[] values;
    delete[] indices;
    delete[] offset;
    delete[] ptr;
}


template<typename T>
void SparseMatrix<T>::clear(){
    fill(values, values + nnz, T());
}

template<typename T>
void SparseMatrix<T>::printSparse(const unsigned &spacing) const {
    for(unsigned i = 0; i < height; i++){
        for(unsigned j = 0; j < offset[i + 1] - offset[i]; j++)
            cout << setw(spacing) << values[offset[i] + j] << " (" << i << ", " << indices[offset[i] + j] << ")" << endl; 
    }
}

template<typename T>
void SparseMatrix<T>::printDense(const unsigned &spacing, const unsigned &precision) const {
    unsigned position;
    for(unsigned i = 0; i < height; i++){
        position = 0;
        for(unsigned j = 0; j < offset[i + 1] - offset[i]; j++){
            while(position != indices[offset[i] + j]){
                cout << setw(spacing) << 0;
                position++;
            }
            cout << setw(spacing) << setprecision(precision) << values[offset[i] + j];
            position++;
        }
        for(;position < width; position++){
            cout << setw(spacing) << 0;
        }
        cout << "\n";
    }
}
template<typename T>
vector<vector<Node<T>>> SparseMatrix<T>::getData() const {
    vector<vector<Node<T>>> data(height);
    for(unsigned i = 0; i < height; i++){
        for(unsigned j = 0; j < offset[i + 1] - offset[i]; j++){
            data[i].push_back({values[offset[i] + j], indices[offset[i] + j]});
        }
    }
    return data;
}


template<typename T>
SparseMatrix<T> SparseMatrix<T>::getTransposed() const {
    SparseMatrix<T> transposed;
    transposed.height = width;
    transposed.width = height;
    transposed.nnz = nnz;
    transposed.values = new T[nnz];
    transposed.indices = new unsigned[nnz];
    transposed.ptr = new unsigned[nnz];
    transposed.offset = new  unsigned[width + 1];

    vector<unsigned> length_sizes(width, 0);
    for(unsigned i = 0; i < nnz; i++){
        length_sizes[indices[i]]++;
    }
    transposed.offset[0] = 0;
    for(unsigned i = 0; i < width; i++){
        transposed.offset[i + 1] = transposed.offset[i] + length_sizes[i];
    }
    fill(length_sizes.begin(), length_sizes.end(), 0);

    for(unsigned i = 0; i < height; i++){
        for(unsigned j = 0; j < offset[i + 1] - offset[i]; j++){
            unsigned index = indices[offset[i] + j];
            transposed.values[transposed.offset[index] + length_sizes[index]] = values[offset[i] + j];
            transposed.indices[transposed.offset[index] + length_sizes[index]] = i;            
            length_sizes[index]++;
        }
    }
    return transposed;
}
template<typename T>
SparseMatrix<T> SparseMatrix<T>::getSquared() const {
    SparseMatrix<T> squared(*this);
    for(unsigned i = 0; i < nnz; i++){
        squared.values[i] = squared.values[i] * squared.values[i];
    }
    return squared;
}
template<typename T>
SparseMatrix<T> SparseMatrix<T>::getInverse() const {
    SparseMatrix<T> inverse(*this);
    for(unsigned i = 0; i < nnz; i++){
        inverse.values[i] = 1.0/inverse.values[i];
    }
    return inverse;
}

template class SparseMatrix<float>;
template class SparseMatrix<double>;