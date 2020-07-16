#ifndef ANALYTICS_H
#define ANALYTICS_H

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/SPQRSupport>
#include "SparseMatrix.hpp"
#include <algorithm>

class Analytics {
    Eigen::SparseMatrix<double> A;
public:
    Analytics(const SparseMatrix<double> A): A(A.getHeight(), A.getWidth())
    {
        Analytics::A.reserve(A.getNNZ());
        vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(A.getNNZ());
        // Populate H
        for(unsigned i = 0; i < A.getHeight(); i++){
            for(unsigned j = 0; j < A.getRowSize(i); j++){
                triplets.push_back(Eigen::Triplet<double>(i, A.getIndex(i, j), A(i, j)));
            }
        }
        Analytics::A.setFromTriplets(triplets.begin(), triplets.end());
    
    }

    int getRank() const {
        Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        return solver.rank(); 
    }

    template<typename T>
    static double computeWrss(vector<T> marginals, const SparseMatrix<T> &H, vector<T> means, vector<T> variances){
        double wrss = 0.0;

        // #pragma omp parallel for reduction(+:wrss)
        for(unsigned i = 0; i < H.getHeight(); i++){
            double sum = 0.0;
            for(unsigned j = 0; j < H.getRowSize(i); j++){
                sum += (H(i,j) * marginals[H.getIndex(i, j)]);
            }
            wrss += (sum - means[i])*(sum - means[i])/variances[i];
        }
        return wrss;
    }

    public:
    template<typename T>
    static double spectralRadius(SparseMatrix<T> &H){
        Eigen::SparseMatrix<double> H_eigen(H.getHeight(), H.getWidth());
        H_eigen.reserve(H.getNNZ());
        vector<Eigen::Triplet<T>> triplets;
        triplets.reserve(H.getNNZ());

        for(unsigned i = 0; i < H.getHeight(); i++){
            for(unsigned j = 0; j < H.getRowSize(i); j++){
                triplets.push_back(Eigen::Triplet<T>(i, H.getIndex(i, j), H(i, j)));
            }
        }
        H_eigen.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::EigenSolver<Eigen::SparseMatrix<double>> eigen_solver(H_eigen);
        Eigen::VectorXd eigen_values = eigen_solver.eigenvalues(); 
        unsigned max_index = 0;
        for(unsigned i = 0; i < eigen_values.size(); i++){
            if(abs(eigen_values(i)) > abs(eigen_values(i))){
                max_index = i;
            }
        }
        return abs(eigen_values(max_index));
    }
};

#endif // ANALYTICS_H