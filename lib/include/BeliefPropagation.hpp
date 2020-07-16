#ifndef BELIEF_PROPAGATION_H
#define BELIEF_PROPAGATION_H

#include "SparseMatrix.hpp"
#include <vector>

template<typename T>
class BeliefPropagation {
    static constexpr double damphing_probability = 0.6;
    static constexpr double alpha = 0.4;
    static constexpr unsigned dampless_count = 10;
    bool damphing_enabled = false;

    SparseMatrix<T> H;
    vector<T> means, variances, direct_means, direct_variances;
    SparseMatrix<T> m_variable_to_factor, m_factor_to_variable, v_variable_to_factor, v_factor_to_variable;
    unsigned* damphing;
    
public:
    BeliefPropagation(SparseMatrix<T> &H, vector<T> &means, vector<T> &variances, vector<T> &direct_means, vector<T> &direct_variances);
    void initilizeMessages();
    void initilizeDamping();
    void linkMessages();
    
    void sendVariableToFactor(const unsigned &it);
    void sendVariableToFactor_mean(const unsigned &it);

    void sendFactorToVariable(const unsigned &it);
    void sendFactorToVariable_mean(const unsigned &it);

    void run(const unsigned &iteration_count = 200, const unsigned &iteration_count_variance = 20);
    const vector<T> computeMarginals() const;

    //Getters
    SparseMatrix<T>& getMeanVariableToFactor() { return m_variable_to_factor; }
    SparseMatrix<T>& getVarianceVariableToFactor() { return v_variable_to_factor; }
    SparseMatrix<T>& getMeanFactorToVariable() { return m_factor_to_variable; }
    SparseMatrix<T>& getVarianceFactorToVariable() { return v_factor_to_variable; }
    const unsigned* getDamphing() const { return damphing; }

    ~BeliefPropagation(){ delete[] damphing; }
};
#endif // BELIEF_PROPAGATION_H