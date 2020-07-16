#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <vector>
#include <iostream>
using namespace std;

class Equations {

public:
    // Calculate sum
    template<typename T>
    static const T sumMeanVariableToFactor(const T* vector_v_factor_to_variable, const T* vector_m_factor_to_variable, const unsigned size){
        T sum = T();
        for(unsigned i = 0; i < size; i++){
            sum +=  vector_m_factor_to_variable[i]/vector_v_factor_to_variable[i];
        }
        return sum;
    }

    template<typename T>
    static const T sumVarianceVariableToFactor(const T* vector_v_factor_to_variable, const unsigned size){
        T sum = T();
        for(unsigned i = 0; i < size; i++){
            sum += 1.0/vector_v_factor_to_variable[i];
        }
        return sum;
    }
   
    template<typename T>
    static const T sumMeanFactorToVariable(const T* vector_factor, const T* vector_m_variable_to_factor, const unsigned size){
        T sum = T();
        for(unsigned i = 0; i < size; i++){
            sum += vector_factor[i]*vector_m_variable_to_factor[i];
        }
        return sum;
    }
    template<typename T>
    static const T sumVarianceFactorToVariable(const T* vector_factor, const T* vector_v_variable_to_factor, const unsigned size){
        T sum = T();
        for(unsigned i = 0; i < size; i++){
            const T factorSquared = vector_factor[i]*vector_factor[i];
            sum += factorSquared*vector_v_variable_to_factor[i];
        }
        return sum;
    }

    // Kahan correction sum
    template<typename T>
    static const T kahan_sumMeanVariableToFactor(const T* vector_v_factor_to_variable, const T* vector_m_factor_to_variable, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = vector_m_factor_to_variable[i]/vector_v_factor_to_variable[i] - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t; 
        }
        return sum;
    }
    template<typename T>
    static const T kahan_sumVarianceVariableToFactor(const T* vector_v_factor_to_variable, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = 1.0/vector_v_factor_to_variable[i] - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t; 
        }
        return sum;
    }
    template<typename T>
    static const T kahan_sumMeanFactorToVariable(const T* vector_factor, const T* vector_m_variable_to_factor, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = vector_factor[i]*vector_m_variable_to_factor[i] - c;
            
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t;
        }
        return sum;
    }
    template<typename T>
    static const T kahan_sumVarianceFactorToVariable(const T* vector_factor, const T* vector_v_variable_to_factor, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T factorSquared = vector_factor[i]*vector_factor[i];
            const T next = factorSquared*vector_v_variable_to_factor[i] - c;
            const T t = sum + next;
            c = (t - sum) - next;
            sum = t;
        }
        return sum;
    }

    // Neumaier correction sum
    template<typename T>
    static const T neumaier_sumMeanVariableToFactor(const T* vector_v_factor_to_variable, const T* vector_m_factor_to_variable, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = vector_m_factor_to_variable[i]/vector_v_factor_to_variable[i];
            const T t = sum + next;
            if(abs(sum) >= abs(next)){
                c += (sum - t) + next;
            } else {
                c += (next - t) + sum; 
            }
            sum = t; 
        }
        return sum + c;
    }
    template<typename T>
    static const T neumaier_sumVarianceVariableToFactor(const T* vector_v_factor_to_variable, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = 1.0/vector_v_factor_to_variable[i];
            const T t = sum + next;
            if(abs(sum) >= abs(next)){
                c += (sum - t) + next;
            } else {
                c += (next - t) + sum;
            }
            sum = t; 
        }
        return sum + c;
    }
        template<typename T>
    static const T neumaier_sumMeanFactorToVariable(const T* vector_factor, const T* vector_m_variable_to_factor, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T next = vector_factor[i]*vector_m_variable_to_factor[i];
            const T t = sum + next;
            if(abs(sum) >= abs(next)){
                c += (sum - t) + next;
            } else {
                c += (next - t) + sum;
            }
            sum = t;
        }
        return sum + c;
    }
    template<typename T>
    static const T neumaier_sumVarianceFactorToVariable(const T* vector_factor, const T* vector_v_variable_to_factor, const unsigned size){
        T sum = T();
        T c = T();
        for(unsigned i = 0; i < size; i++){
            const T factorSquared = vector_factor[i]*vector_factor[i];
            const T next = factorSquared*vector_v_variable_to_factor[i];
            const T t = sum + next;
            if(abs(sum) >= abs(next)){
                c += (sum - t) + next;
            } else {
                c += (next - t) + sum;
            }
            sum = t;
        }
        return sum + c;
    }
    // Calculate individual messages
    template<typename T>
    static inline const T calculateMeanVariableToFactor(const T& sum, const T& v_variable_to_factor, const T& v_factor_to_variable, const T& m_factor_to_variable, const T& m_direct_factor_to_variable){
        return v_variable_to_factor*(sum - m_factor_to_variable/v_factor_to_variable + m_direct_factor_to_variable);
    }   
    template<typename T>
    static inline const T calculateVarianceVariableToFactor(const T& sum, const T& v_factor_to_variable,  const T& v_direct_factor_to_variable){
        return 1.0/(sum - 1.0/v_factor_to_variable + v_direct_factor_to_variable);
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
