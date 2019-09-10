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
SparseMatrix<T>::SparseMatrix(string file_name) {
    std::ifstream fileInput(file_name);
    if(fileInput.fail()){
        throw domain_error("Failed opening file in Sparse Matrix constructor");
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
}
template<typename T>
SparseMatrix<T>::SparseMatrix(SparseMatrix &sm):
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


template class SparseMatrix<float>;
template class SparseMatrix<double>;