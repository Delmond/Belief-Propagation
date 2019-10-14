#include <iostream>
#include <string>
#include <chrono>

#include "lib/include/SparseMatrix.hpp"
#include "lib/include/VectorIO.hpp"
#include "lib/include/BeliefPropagation.hpp"
#include "lib/include/Preprocessing.hpp"
#include "lib/include/Analytics.hpp"
#include "lib/include/WLS.hpp"

using namespace std;


int main(int argc, char* argv[]){
    string folder = (argc == 1)? "toy": argv[1];

    string H_path = "bin/" + folder + "/factors.txt";
    string means_path = "bin/" + folder + "/means.txt";
    string variances_path = "bin/" + folder + "/variances.txt";

    SparseMatrix<double> H = SparseMatrix<double>::fromFileCSC("bin/data283803x70000/factors_csc.txt");
    vector<double> means = VectorIO::loadFromFile<double>("bin/data283803x70000/means.txt");
    vector<double> variances = VectorIO::loadFromFile<double>("bin/data283803x70000/variances.txt");

    // SparseMatrix<double> H = SparseMatrix<double>::fromFileDense(H_path);
    // vector<double> means = VectorIO::loadFromFile<double>(means_path);
    // vector<double> variances = VectorIO::loadFromFile<double>(variances_path);
    cout << "Preprocessing - Reordering nodes" << endl;

    Preprocessing::reorderNodes<double>(H, means, variances);
    cout << "Preprocessing - Adding Virtual nodes" << endl;

    Preprocessing::addVirtualNodes<double>(H, means, variances);
    cout << "Preprocessing - Normalizing variances" << endl;
    Preprocessing::normalizeVariances<double>(variances);
    cout << "Preprocessing - done" << endl;

    cout<< "H: number of factor nodes: " << H.getHeight() << ", number of varaible nodes: " << H.getWidth() << endl;
    // H.getTransposed();

    BeliefPropagation<double> beliefPropagation(H, means, variances);     


    beliefPropagation.initilizeMessages();
    beliefPropagation.initilizeDamping();
    beliefPropagation.linkMessages();

    // cout << "Initial means - variable to factor" << endl;
    // beliefPropagation.getMeanVariableToFactor().printDense(6);
    // cout << "Initial variances - variable to factor" << endl;
    // beliefPropagation.getVarianceVariableToFactor().printDense(6);
    
    // cout << "Initial means - factor to variable" << endl;
    // beliefPropagation.getMeanFactorToVariable().printDense(6);
    // cout << "Initial variances - factor to variable" << endl;
    // beliefPropagation.getVarianceFactorToVariable().printDense(6);

    auto start = chrono::high_resolution_clock::now();
    beliefPropagation.run(200);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Iteration time: " << duration.count() << "s" << endl;
    vector<double> marginals = beliefPropagation.computeMarginals();

    // for(int i = 0; i < 2000; i++){
    //     beliefPropagation.run(1);
    //     vector<double> marginals = beliefPropagation.computeMarginals();
    //     double wrss = Analytics::computeWrss<double>(marginals, H, means, variances);
    //     cout << wrss << endl;
    // }
    // cout << "Marginals: " << endl;
    // for(unsigned i = 0; i < marginals.size(); i++){
    //     cout << marginals[i] << " ";
    // }
    // cout << endl;
    // WLS<double> wls(H, means, variances);
    // wls.run_QR();
    // vector<double> wls_marginals = wls.computeMarginals();

    // cout << "WLS Marginals: " << endl;
    // for(unsigned i = 0; i < wls_marginals.size(); i++){
    //     cout << wls_marginals[i] << " ";
    // }
    cout << endl;
    return 0;
}