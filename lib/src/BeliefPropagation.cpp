#include "../include/BeliefPropagation.hpp"
#include "../include/Equations.hpp"
#include "omp.h"
#include <iostream>
#include <random>

using namespace std;

template<typename T>
BeliefPropagation<T>::BeliefPropagation(SparseMatrix<T> H, vector<T> means, vector<T> variances):
    H(H),
    H_inverse(H.getInverse()),
    H_squared(H.getSquared()),
    H_inverse_squared(H.getSquared().getInverse()),
    H_transposed(H.getTransposed()),
    means(means),
    variances(variances),
    m_variable_to_factor(H),
    m_factor_to_variable(H.getTransposed()),
    v_variable_to_factor(H),
    v_factor_to_variable(H.getTransposed())
    {}

template<typename T>
void BeliefPropagation<T>::initilizeMessages(){
    m_factor_to_variable.clear();
    v_factor_to_variable.clear();
    for(unsigned i = 0; i < H.getWidth(); i++) {
        m_factor_to_variable(i, 0) = means[i];
        v_factor_to_variable(i, 0) = variances[i];
    }
    m_variable_to_factor.clear();
    v_variable_to_factor.clear();
    for(unsigned i = H.getWidth(); i < H.getHeight(); i++){
        for(unsigned j = 0 ; j < m_variable_to_factor.getRowSize(i); j++){
            m_variable_to_factor(i, j) = means[j];
            v_variable_to_factor(i, j) = 1.0/variances[j];
        }
    }
}
template<typename T>
void BeliefPropagation<T>::initilizeDamping() {
    default_random_engine generator;
    bernoulli_distribution distribution(damphing_probability);
    for(unsigned i = 0; i < H.getHeight(); i++){
        damphing.push_back(vector<int>(H.getRowSize(i), 0));
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            damphing[i][j] = (distribution(generator))? 1 : 0;
        }
    }
    damphing_enabled = true;
}


template<typename T>
void BeliefPropagation<T>::linkMessages(){
    vector<unsigned> indices(m_factor_to_variable.getHeight(), 0);
    unsigned index = 0;
    for(unsigned i = 0; i < m_variable_to_factor.getHeight(); i++){
        for(unsigned j = 0; j < m_variable_to_factor.getRowSize(i); j++){
            index = m_variable_to_factor.getIndex(i, j);
            m_factor_to_variable.link(index, indices[index], m_variable_to_factor.getAddress(i, j));
            v_factor_to_variable.link(index, indices[index], v_variable_to_factor.getAddress(i, j));
            
            m_variable_to_factor.link(i, j, m_factor_to_variable.getAddress(index, indices[index]));
            v_variable_to_factor.link(i, j, v_factor_to_variable.getAddress(index, indices[index]));
            
            indices[index]++;
        }
    }
}


template<typename T>
void BeliefPropagation<T>::sendVariableToFactor(){
    #pragma omp parallel for
    for(unsigned i = H.getWidth(); i <  H.getHeight(); i++){
        T mean_sum = Equations::correctionsumMeanFactorToVariable<T>(H.getVector(i), m_variable_to_factor.getVector(i));
        T variance_sum = Equations::correctionsumVarianceFactorToVariable<T>(H.getVector(i), v_variable_to_factor.getVector(i));
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            const T computedVariance = Equations::calculateVarianceFactorToVariable(variance_sum, H(i, j), variances[i], v_variable_to_factor(i, j));
            const T computedMean = Equations::calculateMeanFactorToVariable(mean_sum, H(i, j), means[i], m_variable_to_factor(i, j));
            cout << "mean_sum: " << mean_sum << ", variance_sum: " << variance_sum << endl;
       
            *v_variable_to_factor.getLink(i, j) = computedVariance;
            // if(damphing_enabled){
                *m_variable_to_factor.getLink(i, j) = (1 - damphing[i][j])*computedMean + damphing[i][j]*(alpha*computedMean + (*m_variable_to_factor.getLink(i, j)) * (1 - alpha));
            // } else {
            //     *m_variable_to_factor.getLink(i, j) = computedMean;
            // }
            
        }
    }
}

template<typename T>
void BeliefPropagation<T>::sendFactorToVariable(){
    #pragma omp parallel for
    for(unsigned i = 0; i < H_transposed.getHeight(); i++){
        T mean_sum = Equations::correctionsumMeanVariableToFactor(v_factor_to_variable.getVector(i), m_factor_to_variable.getVector(i));
        T variance_sum = Equations::correctionsumVarianceVariableToFactor(v_factor_to_variable.getVector(i));
        cout << "mean_sum: " << mean_sum << ", variance_sum: " << variance_sum << endl;
       
        for(unsigned j = 1; j < H_transposed.getRowSize(i); j++){
            const T computedVariance = Equations::calculateVarianceVariableToFactor(variance_sum, v_factor_to_variable(i, j));
            const T computedMean = Equations::calculateMeanVariableToFactor(mean_sum, computedVariance, v_factor_to_variable(i, j), m_factor_to_variable(i, j));
            *v_factor_to_variable.getLink(i, j) = computedVariance;
            *m_factor_to_variable.getLink(i, j) = computedMean;
        }
    }
}

template<typename T>
void BeliefPropagation<T>::run(const unsigned &iteration_count){
   for(unsigned i = 0; i < iteration_count; i++){
        cout << "Iteration " << i << endl;
        cout << "varinace -> factor " << endl; 
        sendVariableToFactor();
        cout << "factor -> varinace " << endl; 
        sendFactorToVariable();
    }
}
template<typename T>
const vector<T> BeliefPropagation<T>::computeMarginals() const {
    vector<T> marginals(H.getWidth(), T());
    vector<T> variances_final(H.getWidth(), T());
    for(unsigned i = 0; i < H.getWidth(); i++){
        for(unsigned j = 0; j < v_factor_to_variable.getRowSize(i); j++){
            variances_final[i] += 1.0/v_factor_to_variable(i, j);
        }
        variances_final[i] = 1.0/variances_final[i];
    }
    
    for(unsigned i = 0; i < H.getWidth(); i++){
        for(unsigned j = 0; j < m_factor_to_variable.getRowSize(i); j++){
            marginals[i] += m_factor_to_variable(i, j)/v_factor_to_variable(i, j); 
        }
        marginals[i] *= variances_final[i];
  }
  return marginals;
}


template class BeliefPropagation<float>;
template class BeliefPropagation<double>;
