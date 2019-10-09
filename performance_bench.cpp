#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>
#include "lib/include/SparseMatrix.hpp"
#include "lib/include/VectorIO.hpp"
#include "lib/include/BeliefPropagation.hpp"
#include "lib/include/Preprocessing.hpp"



double bench_iteration(SparseMatrix<double> &H, vector<double> &means, vector<double> &variances, int iteration_count){
    BeliefPropagation<double> beliefPropagation(H, means, variances);     

    beliefPropagation.initilizeMessages();
    beliefPropagation.initilizeDamping();
    beliefPropagation.linkMessages();
    
    auto start = chrono::high_resolution_clock::now();
    beliefPropagation.run(iteration_count);
    auto end = chrono::high_resolution_clock::now();
    
    chrono::duration<double> duration = end - start;
    return duration.count();
}


int main() {
    constexpr int initial_iteration_count = 5000;
    constexpr int iteration_increment = 5000;
    constexpr int repeat = 10;

    string H_path = "bin/damp/factors.txt";
    string means_path = "bin/damp/means.txt";
    string variances_path = "bin/damp/variances.txt";

    cout << "Bench testing started" << endl;
    cout << "Loading damphing dataset" << endl;
    
    SparseMatrix<double> H(H_path);
    vector<double> means = VectorIO::loadFromFile<double>(means_path);
    vector<double> variances = VectorIO::loadFromFile<double>(variances_path);
    
    cout << "Loaded damphing dataset H " << H.getHeight() << "x" << H.getWidth() << endl;

    cout << "Preprocessing started..." << endl;

    Preprocessing::reorderNodes<double>(H, means, variances);
    Preprocessing::addVirtualNodes<double>(H, means, variances);
    Preprocessing::normalizeVariances<double>(variances);
    
    cout << "Preprocessing done. New datase H size " << H.getHeight() << "x" << H.getWidth() << endl;
    vector<vector<double>> results(10, vector<double>());

    cout << "Bench settings:" << endl;
    cout << "\tinitial iteration count: " << initial_iteration_count << endl;
    cout << "\titeration increment: " << iteration_increment << endl;
    cout << "\tall benchmarks will repeat: " << repeat << " times" << endl;
    for(int t = 0; t < repeat; t++){
        cout << "Test repeation: " << t + 1 << endl;
        for(int i = 0; i < 10; i++){
            cout << "\tTesting for " << initial_iteration_count + i*iteration_increment << " iterations" << endl;
            double elapsed_time = bench_iteration(H, means, variances, initial_iteration_count + i*iteration_increment);
            cout << "Time: " << elapsed_time << "s" << endl;
            results[i].push_back(elapsed_time);
        }
    }

    string output_name = "4cores.txt";
    ofstream output(output_name);

    for(int i = 0; i < 10; i++){
        double result_mean = 0;
        for(int j = 0; j < repeat; j++){
            result_mean += results[i][j];
        }
        result_mean /= repeat;
        output << initial_iteration_count + i*iteration_increment << ", " << setprecision(15) << result_mean << endl;
    }

    return 0;
}