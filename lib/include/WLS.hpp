#ifndef WLS_H
#define WLS_H

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/SPQRSupport>
// #include <eigen3/Eigen/PardisoSupport>
#include <eigen3/Eigen/CholmodSupport>
#include <eigen3/Eigen/SuperLUSupport>
#include <eigen3/Eigen/UmfPackSupport>

#include "SparseMatrix.hpp"

template<typename T>
class WLS {
    Eigen::SparseMatrix<T> H, means, W, W_sqrt;
    Eigen::Matrix<T, Eigen::Dynamic, 1> marginals;
  
public:
    

    WLS(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances);
    WLS(Eigen::SparseMatrix<T> &H, vector<T> &means, vector<T> &variances);
    void run_mldivide();
    void run_QR();
    void run_SPQR();
    void run_LLT();
    void run_LU();
    void run_LDLT();
    void run_BiCGSTAB();
    void run_ConjugateGradient();
    void run_superLU();
    void run_umfpackLU();
    void print_marginals(unsigned start = 0, unsigned end = 0);
    vector<T> computeMarginals();
};

#endif // WLS_H