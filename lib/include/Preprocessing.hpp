#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "SparseMatrix.hpp"
#include <algorithm>

class Preprocessing {

public:

    template<typename T>
    static void reorderNodes(SparseMatrix<T> &factors, vector<T> &means, vector<T> &variances){
        unsigned i = 0;
        while(i != factors.getWidth()){
            unsigned next_index = (factors.getRowSize(i) == 1)? factors.getIndex(i, 0) : factors.getHeight();
            for(unsigned j = i + 1; j < factors.getHeight(); j++){
                if(factors.getRowSize(j) != 1)
                    continue;
                if(next_index == factors.getHeight() || factors.getIndex(j, 0) < factors.getIndex(next_index, 0)){
                    next_index = j;
                }

            }
            if(next_index == factors.getHeight())
                return;
            if(factors.getRowSize(i) != 1 || factors.getIndex(next_index, 0) < factors.getIndex(i, 0)){
                factors.swapRows(i, next_index);
                swap(means[i], means[next_index]);
                swap(variances[i], variances[next_index]);
            }
            i++;
        }
    }

    template<typename T>
    static void addVirtualNodes(SparseMatrix<T> &factors, vector<T>& means, vector<T>& variances){
        const T virtualMean = T();
        const T virtualVariance = T(1e2);
        unsigned i = 0;
        while(i < factors.getWidth()){
            if(factors.getRowSize(i) == 1 && factors.getIndex(i, 0) == i) {
                i++;
                continue;
            }
            vector<Node<T>> row(1, Node<T>({1, i, nullptr}));
            factors.insertRow(row, i);
            means.insert(means.begin() + i, virtualMean);
            variances.insert(variances.begin() + i, virtualVariance);
        }
    }

    template<typename T>
    static void normalizeVariances(vector<T> &variances){
        for(unsigned i = 0; i < variances.size(); i++){
            if(variances[i] < 1e-8)
                variances[i] = 1e-8;
            if(variances[i] > 1e2)
                variances[i] = 1e2;
        }
    }

};

#endif // PREPROCESSING_H