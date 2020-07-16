#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/CholmodSupport>

#include "lib/include/SparseMatrix.hpp"
#include "lib/include/VectorIO.hpp"
#include "lib/include/BeliefPropagation.hpp"
#include "lib/include/Preprocessing.hpp"
#include "lib/include/Analytics.hpp"
#include "lib/include/WLS.hpp"
#include "lib/include/Benchmarking.hpp"
using namespace std;

#define EIGEN_USE_BLAS

int main(int argc, char* argv[]){
    Eigen::initParallel();

    string folder = (argc == 1)? "data14x12": argv[1];

    string H_path = "bin/" + folder + "/factors.txt";
    string means_path = "bin/" + folder + "/means.txt";
    string variances_path = "bin/" + folder + "/variances.txt";

    // string csv_folder = "data33x14";
    // string csv_folder = "data283803x70000";
    string csv_folder = "data69999x69999";
    // string csv_folder = "data5997x2000";
    // string csv_folder = "sherman2";
    // string csv_folder = "area20";

    SparseMatrix<double> H = SparseMatrix<double>::fromFileCSC("bin/" + csv_folder + "/factors_csc.txt");
    vector<double> means = VectorIO::loadFromFile<double>("bin/" + csv_folder + "/means.txt");
    vector<double> variances = VectorIO::loadFromFile<double>("bin/" + csv_folder + "/variances.txt");

    // SparseMatrix<double> H = SparseMatrix<double>::fromFileDense(H_path);
    // vector<double> means = VectorIO::loadFromFile<double>(means_path);
    // vector<double> variances = VectorIO::loadFromFile<double>(variances_path);


    Analytics analytics(H);
    cout << "Information about A:" << endl;
    cout << "Height: " << H.getHeight() << ", Width: " << H.getWidth() << endl;
    cout << "Rank: " << analytics.getRank() << endl;
    cout << endl;

    SparseMatrix<double> H_bp(H);
    vector<double> means_bp = means, variances_bp = variances, direct_means, direct_variances;

    auto start = chrono::high_resolution_clock::now();
    Preprocessing::removeDirectNodes(H_bp, means_bp, variances_bp, direct_means, direct_variances);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    
    cout << "Preprocessing - done in: " << duration.count() <<  endl;

    BeliefPropagation<double> beliefPropagation(H_bp, means_bp, variances_bp, direct_means, direct_variances);     
    beliefPropagation.initilizeMessages();
    beliefPropagation.initilizeDamping();
    beliefPropagation.linkMessages();    

    start = chrono::high_resolution_clock::now();
    beliefPropagation.run(6000, 6000);
    end = chrono::high_resolution_clock::now();
    duration = end - start;
    vector<double> marginals = beliefPropagation.computeMarginals();

    #ifdef DEBUG
    cout << "Marginals: ";
    for(unsigned i = 0; i < marginals.size(); i++){
        cout << marginals[i] << " ";
    }
    cout << endl;
    #endif
    cout << "BP iteration time: " << duration.count() << "s" << endl;
    
    // for(int i = 0; i < 2000; i++){
    //     beliefPropagation.run(1);
    //     vector<double> marginals = beliefPropagation.computeMarginals();
    //     double wrss = Analytics::computeWrss<double>(marginals, H, means, variances);
    //     cout << wrss << endl;
    // }
    // cout << endl;

    // WLS<double> wls(H, means, variances);
    // wls.run_LLT();
    // wls.computeMarginals();
    // auto wls_marginals = wls.computeMarginals();
    // for(unsigned i = 0; i < wls_marginals.size(); i++){
    //     double diff =  abs(marginals[i] - wls_marginals[i]);
    //     cout << diff << endl;

    // }
    // start = chrono::high_resolution_clock::now();
    // wls.run_ConjugateGradient();
    // // wls.run_QR();
    // end = chrono::high_resolution_clock::now();
    // duration = end - start;
    // cout << "Iteration time: " << duration.count() << "s" << endl;



    // auto QR_bench = Benchmarking::benchQR(H, means, variances);
    // cout << "QR iteration time: " << QR_bench.first.count() << endl;

    // vector<Eigen::Triplet<double>> triplets;
    // triplets.push_back(Eigen::Triplet<double>(0, 0, 4e60));
    // triplets.push_back(Eigen::Triplet<double>(0, 1, -2e60));
    // triplets.push_back(Eigen::Triplet<double>(1, 0, -2e60));
    // triplets.push_back(Eigen::Triplet<double>(1, 1, 1e60));
    // Eigen::SparseMatrix<double> test;
    // test.setFromTriplets(triplets.begin(), triplets.end());
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(test);

    // Eigen::Vector2d b;
    // b.coeffRef(0) = 1;
    // b.coeffRef(1) = 1;
     

    
    


    // 
    auto LLT_bench = Benchmarking::benchLLT(H, means, variances);
    cout << "LLT iteration time: " << LLT_bench.first.count() << endl;

    auto LU_bench = Benchmarking::benchLLT(H, means, variances);
    cout << "LU iteration time: " << LU_bench.first.count() << endl;

    auto LDLT_bench = Benchmarking::benchLDLT(H, means, variances);
    cout << "LDLT iteration time: " << LDLT_bench.first.count() << endl;

    auto spqr_bench = Benchmarking::benchSPQR(H, means, variances);
    cout << "SPQR iteration time: " << spqr_bench.first.count() << endl;

    auto superlu_bench = Benchmarking::benchSuperLU(H, means, variances);
    cout << "SuperLU iteration time: " << superlu_bench.first.count() << endl;

    auto umfpacklu_bench = Benchmarking::benchSuperLU(H, means, variances);
    cout << "UmfpackLU iteration time: " << umfpacklu_bench.first.count() << endl;


    double Wrss_BP = Analytics::computeWrss(marginals, H, means, variances);
    // double Wrss_QR = QR_bench.second;
    double Wrss_SPQR = spqr_bench.second;
    double Wrss_LLT = LLT_bench.second;
    double Wrss_LU = LU_bench.second;
    double Wrss_LDLT = LDLT_bench.second;
    double Wrss_Super_LU = superlu_bench.second;
    double Wrss_Umfpack_LU = umfpacklu_bench.second;
    

    cout << endl;
    cout << "BP WRSS: " << Wrss_BP << endl;
    // cout << "QR WRSS: " << Wrss_QR << endl;
    cout << "SPQR WRSS: " << Wrss_SPQR << endl;
    cout << "LLT WRSS: " << Wrss_LLT << endl;
    cout << "LU WRSS: " << Wrss_LU << endl;
    cout << "LDLT WRSS: " << Wrss_LDLT << endl;
    cout << "SuperLU WRSS: " << Wrss_Super_LU << endl;
    cout << "UmfpackLU WRSS: " << Wrss_Umfpack_LU << endl;
    
    cout << "Ratio: " << Wrss_BP/Wrss_LLT << endl;
    // cout << "Spectral radius: " << Analytics::spectralRadius(H) << endl;



    return 0;
}