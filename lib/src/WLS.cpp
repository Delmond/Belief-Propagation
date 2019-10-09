#include "../include/WLS.hpp"

#include <iostream> 

template<typename T>
WLS<T>::WLS(SparseMatrix<T> H, vector<T> means, vector<T> variances):
    H(H.getHeight(), H.getWidth()),
    means(means.size(), 1),
    W(variances.size(), variances.size())
{
    // Populate H
    for(unsigned i = 0; i < H.getHeight(); i++){
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            WLS::H.insert(i, H.getIndex(i, j)) = H(i, j);
        }
    }
    // Populate W
    for(unsigned i = 0; i < variances.size(); i++){
        WLS::W.insert(i, i) = 1/variances[i];
    }
    // Populate means
    for(unsigned i = 0; i < means.size(); i++){
        WLS::means.insert(i, 0) = means[i];
    }
}
template<typename T>
void WLS<T>::run_mldivide() {
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::SparseMatrix<double> B = H.transpose() * W * means;
}

template<typename T>
void WLS<T>::run_QR(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::SparseMatrix<double> B = H.transpose() * W * means;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);

}
template<typename T>
vector<double> WLS<T>::computeMarginals(){
    vector<double> result(H.cols());
    for(unsigned i = 0; i < result.size(); i++){
        result[i] = marginals.coeff(i, 0);
    }
    return result;
}


template class WLS<float>;
template class WLS<double>;