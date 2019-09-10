#ifndef BELIEF_PROPAGATION_H
#define BELIEF_PROPAGATION_H

#include "SparseMatrix.hpp"
#include <vector>

template<typename T>
class BeliefPropagation {
    SparseMatrix<T> factors;
    vector<T> means, variances;
    SparseMatrix<T> m_state_to_factor, m_factor_to_state, v_state_to_factor, v_factor_to_state;

public:
    BeliefPropagation(SparseMatrix<T> factors, vector<T> means, vector<T> variances);
    void initilizeMessages();
    void linkMessages();
    void sendStateToFactor();
    void sendFactorToState();

    //Getters
    SparseMatrix<T>& getMeanStateToFactor() { return m_state_to_factor; }
    SparseMatrix<T>& getVarianceStateToFactor() { return v_state_to_factor; }
    SparseMatrix<T>& getMeanFactorToState() { return m_factor_to_state; }
    SparseMatrix<T>& getVarianceFactorToState() { return v_factor_to_state; }

};
#endif // BELIEF_PROPAGATION_H