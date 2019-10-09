#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <vector>
#include <iostream>
using namespace std;

class Equations {

public:
    // Calculate sum
    template<typename T>
    static inline const T sumMeanVariableToFactor(const vector<Node<T>> &vector_v_factor_to_variable, const vector<Node<T>> &vector_m_factor_to_variable){
        T sum = T();
        for(unsigned i = 0; i < vector_v_factor_to_variable.size(); i++){
            sum +=  vector_m_factor_to_variable[i].value/vector_v_factor_to_variable[i].value;
        }
        return sum;
    }

    template<typename T>
    static inline const T sumVarianceVariableToFactor(const vector<Node<T>> &vector_v_factor_to_variable){
        T sum = T();
        for(unsigned i = 0; i < vector_v_factor_to_variable.size(); i++){
            sum += 1.0/vector_v_factor_to_variable[i].value;
        }
        return sum;
    }
   
    template<typename T>
    static inline const T sumMeanFactorToVariable(const vector<Node<T>>& vector_factor, const vector<Node<T>> &vector_m_variable_to_factor){
        T sum = T();
        for(unsigned i = 0; i < vector_m_variable_to_factor.size(); i++){
            sum += vector_factor[i].value*vector_m_variable_to_factor[i].value;
        }
        return sum;
    }
    template<typename T>
    static inline const T sumVarianceFactorToVariable(const vector<Node<T>>& vector_factor, const vector<Node<T>> &vector_v_variable_to_factor){
        T sum = T();
        for(unsigned i = 0; i < vector_v_variable_to_factor.size(); i++){
            const T factorSquared = vector_factor[i].value*vector_factor[i].value;
            sum += factorSquared*vector_v_variable_to_factor[i].value;
        }
        return sum;
    }

    // Correction sum
    template<typename T>
    static const T correctionsumMeanVariableToFactor(const vector<Node<T>> &vector_v_factor_to_variable, const vector<Node<T>> &vector_m_factor_to_variable){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < vector_v_factor_to_variable.size(); i++){
            const T next = vector_m_factor_to_variable[i].value/vector_v_factor_to_variable[i].value - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t; 
        }
        return sum;
    }
    template<typename T>
    static const T correctionsumVarianceVariableToFactor(const vector<Node<T>> &vector_v_factor_to_variable){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < vector_v_factor_to_variable.size(); i++){
            const T next = 1.0/vector_v_factor_to_variable[i].value - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t; 
        }
        return sum;
    }
    template<typename T>
    static const T correctionsumMeanFactorToVariable(const vector<Node<T>>& vector_factor, const vector<Node<T>> &vector_m_variable_to_factor){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < vector_m_variable_to_factor.size(); i++){
            const T next = vector_factor[i].value*vector_m_variable_to_factor[i].value - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t;
        }
        return sum;
    }
    template<typename T>
    static const T correctionsumVarianceFactorToVariable(const vector<Node<T>>& vector_factor, const vector<Node<T>> &vector_v_variable_to_factor){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < vector_v_variable_to_factor.size(); i++){
            const T factorSquared = vector_factor[i].value*vector_factor[i].value;
            const T next = factorSquared*vector_v_variable_to_factor[i].value - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t;
        }
        return sum;
    }
    // Calculate individual messages
    template<typename T>
    static inline const T calculateMeanVariableToFactor(const T& sum, const T& v_variable_to_factor, const T& v_factor_to_variable, const T& m_factor_to_variable){
        return v_variable_to_factor*(sum - m_factor_to_variable/v_factor_to_variable);
    }   
    template<typename T>
    static inline const T calculateVarianceVariableToFactor(const T& sum, const T& v_factor_to_variable){
        return 1.0/(sum - 1.0/v_factor_to_variable);
    }

    template<typename T>
    static inline const T calculateMeanFactorToVariable(const T& sum, const T& factor, const T& mean, const T&  m_variable_to_factor){
        const T factorInverse = 1.0/factor;
        return factorInverse * (mean - sum + factor*m_variable_to_factor);
    }

    template<typename T>
    static inline const T calculateVarianceFactorToVariable(const T& sum, const T& factor, const T& variance, const T& v_variable_to_factor){
        const T factorSquared = factor * factor;
        const T factorSquaredInverse = 1.0/factorSquared;
       return factorSquaredInverse * (variance + (sum - factorSquared*v_variable_to_factor));
    }

};


#endif
