#include <iostream>
#include <string>

#include "lib/include/SparseMatrix.hpp"
#include "lib/include/VectorIO.hpp"
#include "lib/include/BeliefPropagation.hpp"

using namespace std;

int main(int argc, char* argv[]){
    string folder = (argc == 1)? "toy": argv[1];

    string factors_path = "bin/" + folder + "/factors.txt";
    string means_path = "bin/" + folder + "/means.txt";
    string variances_path = "bin/" + folder + "/variances.txt";

    SparseMatrix<double> factors(factors_path);
    vector<double> means = VectorIO::loadFromFile<double>(means_path);
    vector<double> variances = VectorIO::loadFromFile<double>(variances_path);

    // factors.getTransposed().printDense();
    BeliefPropagation<double> beliefPropagation(factors, means, variances);    
    
    beliefPropagation.initilizeMessages();
    beliefPropagation.linkMessages();
    
    beliefPropagation.getMeanStateToFactor().printDense(6);
    beliefPropagation.getMeanFactorToState().printDense(6);
    

    return 0;
}