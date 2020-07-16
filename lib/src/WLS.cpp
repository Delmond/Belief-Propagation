#include "../include/WLS.hpp"
#include <iostream> 
#include <cmath>

#define PRINT_MARGINALS
constexpr unsigned start_index = 0;
constexpr unsigned end_index = 12;

template<typename T>
WLS<T>::WLS(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances):
    H(H.getHeight(), H.getWidth()),
    means(means.size(), 1),
    W(variances.size(), variances.size()),
    W_sqrt(variances.size(), variances.size())
{
    WLS::H.reserve(H.getNNZ());
    vector<Eigen::Triplet<T>> triplets;
    triplets.reserve(H.getNNZ());
    // Populate H
    for(unsigned i = 0; i < H.getHeight(); i++){
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            triplets.push_back(Eigen::Triplet<T>(i, H.getIndex(i, j), H(i, j)));
        }
    }
    WLS::H.setFromTriplets(triplets.begin(), triplets.end());
    // Populate W
    for(unsigned i = 0; i < variances.size(); i++){
        WLS::W.insert(i, i) = 1/variances[i];
        WLS::W_sqrt.insert(i, i) = sqrt(1/variances[i]);
    }


    // Populate means
    for(unsigned i = 0; i < means.size(); i++){
        WLS::means.insert(i, 0) = means[i];
    }
}
template<typename T>
WLS<T>::WLS(Eigen::SparseMatrix<T> &H, vector<T> &means, vector<T> &variances): 
    means(means.size(), 1),
    W(variances.size(), variances.size()),
    W_sqrt(variances.size(), variances.size())
{
    WLS::H = H;
    for(unsigned i = 0; i < variances.size(); i++){
        WLS::W.insert(i, i) = 1/variances[i];
        WLS::W_sqrt.insert(i, i) = sqrt(1/variances[i]);
    }
    for(unsigned i = 0; i < means.size(); i++){
        WLS::means.insert(i, 0) = means[i];
    }
}

template<typename T>
void WLS<T>::run_mldivide() {
    Eigen::SparseMatrix<T> A = H.transpose() * W * H;
    Eigen::SparseMatrix<T> B = H.transpose() * W * means;
}

template<typename T>
void WLS<T>::run_SPQR(){
    Eigen::SparseMatrix<double> A = W_sqrt * H;
    Eigen::VectorXd B = W_sqrt * means;
    Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
    // solver.setPivotThreshold(1e-4);
    solver.setSPQROrdering(SPQR_ORDERING_CHOLMOD);
    solver.compute(A);
    WLS::marginals = solver.solve(B);
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}

template<typename T>
void WLS<T>::run_QR(){
    Eigen::SparseMatrix<double> A = W_sqrt * H;
    Eigen::VectorXd B = W_sqrt * means;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
    solver.setPivotThreshold(0);
    solver.compute(A);
    WLS::marginals = solver.solve(B);

}

template<typename T>
void WLS<T>::run_LLT(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}

template<typename T>
void WLS<T>::run_LDLT(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}


template<typename T>
void WLS<T>::run_BiCGSTAB(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
}

template<typename T>
void WLS<T>::run_ConjugateGradient(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
}

template<typename T>
void WLS<T>::run_LU(){
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}
template<typename T>
void WLS<T>::run_superLU(){
    
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
    solver.options().SymmetricMode = yes_no_t::YES;
    solver.compute(A);
    WLS::marginals = solver.solve(B);
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}
template<typename T>
void WLS<T>::run_umfpackLU(){
    
    Eigen::SparseMatrix<double> A = H.transpose() * W * H;
    Eigen::VectorXd B = H.transpose() * W * means;
    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    WLS::marginals = solver.solve(B);    
    #ifdef PRINT_MARGINALS
    print_marginals(start_index, end_index);
    #endif
}

template<typename T>
vector<T> WLS<T>::computeMarginals(){
    vector<T> result(H.cols());
    for(unsigned i = 0; i < result.size(); i++){
        result[i] = marginals.coeff(i, 0);
    }
    return result;
}

template<typename T>
void WLS<T>::print_marginals(unsigned start, unsigned end){
    if(end == 0){
        end = marginals.rows();
    }
    cout << "Marginals: ";
    while(start < end){
        cout << marginals(start) << " ";
        start++;
    }
    cout << endl;
}


// template class WLS<float>;
template class WLS<double>;