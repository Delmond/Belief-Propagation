#ifndef WLS_H
#define WLS_H

#include "../dependancy/Eigen/Sparse"
#include "../dependancy/Eigen/SparseQR"

#include "SparseMatrix.hpp"

template<typename T>
class WLS {
    Eigen::SparseMatrix<double> H, means, W;
    Eigen::SparseMatrix<double> marginals;
  
public:
    WLS(SparseMatrix<T> H, vector<T> means, vector<T> variances);
    void run_mldivide();
    void run_QR();
    void run_LU();
    vector<double> computeMarginals();
};

#endif // WLS_H