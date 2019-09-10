#include "../include/BeliefPropagation.hpp"

template<typename T>
BeliefPropagation<T>::BeliefPropagation(SparseMatrix<T> factors, vector<T> means, vector<T> variances):
    factors(factors),
    means(means),
    variances(variances),
    m_state_to_factor(factors),
    m_factor_to_state(factors.getTransposed()),
    v_state_to_factor(factors),
    v_factor_to_state(factors.getTransposed())
    {}

template<typename T>
void BeliefPropagation<T>::initilizeMessages(){
    for(unsigned i = 0; i < factors.getWidth(); i++) {
        m_state_to_factor(i, 0) = means[i];
        v_state_to_factor(i, 0) = variances[i];
    }
    for(unsigned i = 0; i < factors.getWidth(); i++){
        for(unsigned j = 1 ; j < m_factor_to_state.getRowSize(i); j++){
            m_factor_to_state(i, j) = means[i];
            v_factor_to_state(i, j) = variances[i];
        }
    }

}

template<typename T>
void BeliefPropagation<T>::linkMessages(){
    vector<unsigned> indices(m_factor_to_state.getHeight(), 0);
    unsigned index = 0;
    for(unsigned i = 0; i < m_state_to_factor.getHeight(); i++){
        for(unsigned j = 0; j < m_state_to_factor.getRowSize(i); j++){
            index = m_state_to_factor.getIndex(i, j);
            m_factor_to_state.link(index, indices[index], m_state_to_factor.getAddress(i, j));
            m_state_to_factor.link(i, j, m_factor_to_state.getAddress(index, indices[index]));
            indices[index]++;
        }
    }
}


template<typename T>
void BeliefPropagation<T>::sendStateToFactor(){
    vector<T> sum(m_state_to_factor.getHeight(), T());
    for(unsigned i = 0; i < m_state_to_factor.getHeight(); i++){

    }

}

template<typename T>
void BeliefPropagation<T>::sendFactorToState(){
    for(unsigned i = 0; i < m_factor_to_state.getHeight(); i++){
        
    }
}


template class BeliefPropagation<float>;
template class BeliefPropagation<double>;
