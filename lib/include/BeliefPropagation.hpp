#ifndef BELIEF_PROPAGATION_H
#define BELIEF_PROPAGATION_H

#include "SparseMatrix.hpp"
#include <vector>

template<typename T>
class BeliefPropagation {
    static constexpr double damphing_probability = 0.7;
    static constexpr double alpha = 0.5;
    bool damphing_enabled = false;

    SparseMatrix<T> H, H_inverse, H_squared, H_inverse_squared, H_transposed;
    vector<T> means, variances;
    vector<vector<int>> damphing;
    SparseMatrix<T> m_variable_to_factor, m_factor_to_variable, v_variable_to_factor, v_factor_to_variable;

public:
    BeliefPropagation(SparseMatrix<T> H, vector<T> means, vector<T> variances);
    void initilizeMessages();
    void initilizeDamping();
    void linkMessages();
    void sendVariableToFactor();
    void sendFactorToVariable();
    void run(const unsigned &iteration_count = 200);
    const vector<T> computeMarginals() const;

    //Getters
    SparseMatrix<T>& getMeanVariableToFactor() { return m_variable_to_factor; }
    SparseMatrix<T>& getVarianceVariableToFactor() { return v_variable_to_factor; }
    SparseMatrix<T>& getMeanFactorToVariable() { return m_factor_to_variable; }
    SparseMatrix<T>& getVarianceFactorToVariable() { return v_factor_to_variable; }
    const vector<vector<int>>& getDamphing() const{ return damphing; }

};
#endif // BELIEF_PROPAGATION_H