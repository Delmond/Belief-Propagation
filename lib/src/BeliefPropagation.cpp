#include "../include/BeliefPropagation.hpp"
#include "../include/Equations.hpp"
#include "omp.h"
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

template<typename T>
BeliefPropagation<T>::BeliefPropagation(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances, vector<T> &direct_means, vector<T> &direct_variances):
    H(H),
    means(means),
    variances(variances),
    direct_means(direct_means),
    direct_variances(direct_variances),
    m_variable_to_factor(H),
    m_factor_to_variable(H.getTransposed()),
    v_variable_to_factor(H),
    v_factor_to_variable(H.getTransposed()),
    damphing(nullptr)
    {}

template<typename T>
void BeliefPropagation<T>::initilizeMessages(){
    m_factor_to_variable.clear();
    v_factor_to_variable.clear();
    m_variable_to_factor.clear();
    v_variable_to_factor.clear();

    for(unsigned i = 0; i < H.getHeight(); i++){
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            unsigned index = m_variable_to_factor.getIndex(i, j);
            m_variable_to_factor(i, j) = direct_means[index]/direct_variances[index];
            v_variable_to_factor(i, j) = 1/direct_variances[index];
        }
    }
}
template<typename T>
void BeliefPropagation<T>::initilizeDamping() {
    default_random_engine generator;
    // generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    bernoulli_distribution distribution(damphing_probability);
    damphing = new unsigned[H.getNNZ()];
    for(unsigned i = 0; i < H.getNNZ(); i++){
        damphing[i] = (distribution(generator))? 1 : 0;
    }
    damphing_enabled = true;
}


template<typename T>
void BeliefPropagation<T>::linkMessages(){
    vector<unsigned> length_sizes(m_factor_to_variable.getWidth(), 0);
    for(unsigned i = 0; i < m_factor_to_variable.getHeight(); i++){
        for(unsigned j = 0; j < m_factor_to_variable.getRowSize(i); j++){

            unsigned direct_index = m_factor_to_variable.getDirectIndex(i, j);
            unsigned index = m_factor_to_variable.getIndex(i, j);

            m_variable_to_factor.getLink(index, length_sizes[index]) = direct_index;
            v_variable_to_factor.getLink(index, length_sizes[index]) = direct_index;

            m_factor_to_variable.getLink(i, j) = m_variable_to_factor.getDirectIndex(index, length_sizes[index]);
            v_factor_to_variable.getLink(i, j) = m_variable_to_factor.getDirectIndex(index, length_sizes[index]);
            
            length_sizes[index]++;
        }
    }
    
}


template<typename T>
void BeliefPropagation<T>::sendVariableToFactor(const unsigned &it){


    #pragma omp parallel for
    for(unsigned i = 0; i < H.getHeight(); i++){

        const unsigned size = H.getRowSize(i);

        T mean_sum = Equations::sumMeanFactorToVariable<T>(H.getVector(i), m_variable_to_factor.getVector(i), size);
        T variance_sum = Equations::sumVarianceFactorToVariable<T>(H.getVector(i), v_variable_to_factor.getVector(i), size);

        for(unsigned j = 0; j < H.getRowSize(i); j++){
            const T computedVariance = Equations::calculateVarianceFactorToVariable(variance_sum, H(i, j), variances[i], v_variable_to_factor(i, j));
            const T computedMean = Equations::calculateMeanFactorToVariable(mean_sum, H(i, j), means[i], m_variable_to_factor(i, j));
            unsigned link = v_variable_to_factor.getLink(i, j);
            v_factor_to_variable[link] = computedVariance;

            if(it > dampless_count){
                m_factor_to_variable[link] = (1 - damphing[link])*computedMean + damphing[link]*(alpha*computedMean + (m_factor_to_variable[link]) * (1 - alpha));
            } else {
                m_factor_to_variable[link] = computedMean;
            }    
        }
    }
}
template<typename T>
void BeliefPropagation<T>::sendVariableToFactor_mean(const unsigned &it){
#pragma omp parallel for
    for(unsigned i = 0; i <  H.getHeight(); i++){
        unsigned size = H.getRowSize(i);
        T mean_sum = Equations::sumMeanFactorToVariable<T>(H.getVector(i), m_variable_to_factor.getVector(i), size);
        for(unsigned j = 0; j < H.getRowSize(i); j++){
            const T computedMean = Equations::calculateMeanFactorToVariable(mean_sum, H(i, j), means[i], m_variable_to_factor(i, j));
            unsigned link = m_variable_to_factor.getLink(i, j);
            // if(damphing_enabled){
            if(it > dampless_count){
                m_factor_to_variable[link] = (1 - damphing[link])*computedMean + damphing[link]*(alpha*computedMean + (m_factor_to_variable[link]) * (1 - alpha));
            } else {
                m_factor_to_variable[link] = computedMean;
            }
            
                
        }
    }
}

template<typename T>
void BeliefPropagation<T>::sendFactorToVariable(const unsigned &it){
    #pragma omp parallel for
    for(unsigned i = 0; i < v_factor_to_variable.getHeight(); i++){
        unsigned size = v_factor_to_variable.getRowSize(i);
        T mean_sum = Equations::sumMeanVariableToFactor(v_factor_to_variable.getVector(i), m_factor_to_variable.getVector(i), size);
        T variance_sum = Equations::sumVarianceVariableToFactor(v_factor_to_variable.getVector(i), size);
        for(unsigned j = 0; j < v_factor_to_variable.getRowSize(i); j++){
            const T computedVariance = Equations::calculateVarianceVariableToFactor(variance_sum, v_factor_to_variable(i, j), direct_variances[i]);
            const T computedMean = Equations::calculateMeanVariableToFactor(mean_sum, computedVariance, v_factor_to_variable(i, j), m_factor_to_variable(i, j), direct_means[i]);
            unsigned link = v_factor_to_variable.getLink(i, j);
            v_variable_to_factor[link] = computedVariance;
            m_variable_to_factor[link] = computedMean;
            
            // *v_factor_to_variable.getLink(i, j) = computedVariance;
            // *m_factor_to_variable.getLink(i, j) = computedMean;
        }
    }
}
template<typename T>
void BeliefPropagation<T>::sendFactorToVariable_mean(const unsigned &it){
    #pragma omp parallel for
    for(unsigned i = 0; i < v_factor_to_variable.getHeight(); i++){
        unsigned size = v_factor_to_variable.getRowSize(i);
        T mean_sum = Equations::sumMeanVariableToFactor(v_factor_to_variable.getVector(i), m_factor_to_variable.getVector(i), size);
        for(unsigned j = 0; j < v_factor_to_variable.getRowSize(i); j++){
            unsigned link = m_factor_to_variable.getLink(i, j);
            // const T computedVariance = v_variable_to_factor[link];
            m_variable_to_factor[link] = Equations::calculateMeanVariableToFactor(mean_sum, v_variable_to_factor[link], v_factor_to_variable(i, j), m_factor_to_variable(i, j), direct_means[i]);
            
            // *v_factor_to_variable.getLink(i, j) = computedVariance;
            // *m_factor_to_variable.getLink(i, j) = computedMean;
        }
    }
}


template<typename T>
void BeliefPropagation<T>::run(const unsigned &iteration_count_mean, const unsigned &iteration_count_variance){
    
   for(unsigned i = 0; i < iteration_count_variance; i++){
        sendVariableToFactor(i);
        sendFactorToVariable(i);
    }
    for(unsigned i = iteration_count_variance; i < iteration_count_mean; i++){
        sendVariableToFactor_mean(i);
        sendFactorToVariable_mean(i);    
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
        variances_final[i] = 1.0/(variances_final[i] + direct_variances[i]);
    }
    for(unsigned i = 0; i < H.getWidth(); i++){
        for(unsigned j = 0; j < m_factor_to_variable.getRowSize(i); j++){
            marginals[i] += m_factor_to_variable(i, j)/v_factor_to_variable(i, j);
        }
        marginals[i] =  (marginals[i] + direct_means[i])*variances_final[i];
  }
  return marginals;
}


template class BeliefPropagation<float>;
template class BeliefPropagation<double>;
