#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include<iostream>
#include<chrono>
#include<utility>

#include "./SparseMatrix.hpp"
#include "./WLS.hpp"
#include "./Analytics.hpp"

using namespace std;

class Benchmarking {

public:
    template<typename T>
    static pair<chrono::duration<double>, double> benchQR(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
       
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now(); 
        wls.run_QR();
        auto end = chrono::high_resolution_clock::now();
        
      

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }
    template<typename T>
    static pair<chrono::duration<double>, double> benchSPQR(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_SPQR();
        auto end = chrono::high_resolution_clock::now();

        vector<double> marignals = wls.computeMarginals();

        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }


    template<typename T>
    static pair<chrono::duration<double>, double> benchLLT(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_LLT();
        auto end = chrono::high_resolution_clock::now();
        
        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }

    template<typename T>
    static pair<chrono::duration<double>, double> benchLU(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_LU();
        auto end = chrono::high_resolution_clock::now();

        vector<double> marignals = wls.computeMarginals();
        
        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }

    template<typename T>
    static pair<chrono::duration<double>, double> benchLDLT(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_LDLT();
        auto end = chrono::high_resolution_clock::now();


        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }
    template<typename T>
    static pair<chrono::duration<double>, double> benchSuperLU(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_superLU();
        auto end = chrono::high_resolution_clock::now();
        
        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }

    template<typename T>
    static pair<chrono::duration<double>, double> benchUmfpackLU(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances){
        
        WLS<T> wls(H, means, variances);
        auto start = chrono::high_resolution_clock::now();
        wls.run_umfpackLU();
        auto end = chrono::high_resolution_clock::now();
        
        #ifdef DEBUG
        vector<double> marignals = wls.computeMarginals();
        for(unsigned i = 0; i < marignals.size(); i++){
            cout << marignals[i] << " ";
        }
        cout << endl;
        #endif

        return make_pair(end - start, Analytics::computeWrss<double>(wls.computeMarginals(), H, means, variances));
    }




};

#endif // BENCHMARKING_H