#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <vector>
using namespace std;

class Equations {

public:
    // Calculate sum
    template<typename T>
    static inline const T& sumMeanStateToFactor(const vector<T> &vector_v_factor_to_state, const vector<T> &vector_m_factor_to_state){
        T sum = T();
        for(unsigned i = 0; i < vector_v_factor_to_state.size(); i++){
            sum +=  vector_m_factor_to_state[i]/vector_v_factor_to_state[i];
        }
        return sum;
    }

    template<typename T>
    static inline const T& sumVarianceStateToFactor(const vector<T> &vector_v_factor_to_state){
        T sum = T();
        for(unsigned i = 0; i < vector_v_factor_to_state.size(); i++){
            sum += 1.0/vector_v_factor_to_state[i];
        }
        return sum;
    }

    template<typename T, typename Node_t>
    static inline const T& sumMeanFactorToState(const vector<T>& vector_factor, const vector<Node_t> &vector_m_state_to_factor){
        T sum = T();
        for(unsigned i = 0; i < vector_m_state_to_factor.size(); i++){
            sum += factor[i]*vector_m_state_to_factor[i];
        }
        return sum;
    }
    template<typename T>
    static inline const T& sumVarianceFactorToState(const vector<T>& vector_factor, const vector<T> &vector_v_state_to_factor){
        T sum = T();
        for(unsigned i = 0; i < vector_m_state_to_factor.size(); i++){
            sum += factor[i]*vector_m_state_to_factor[i];
        }
        return sum;
    }

    // Calculate individual messages
    template<typename T>
    static inline const T& calculateMeanStateToFactor(const T& sum, const T& v_factor_to_state, T& m_factor_to_state){
        return sum - m_factor_to_state/v_factor_to_state;
    }   
    template<typename T>
    static inline const T& calculateVarianceStateToFactor(const T& sum, const T& v_factor_to_state, T* output){
        return 1.0/(sum - 1.0/variance);
    }

    template<typename T>
    static inline const T& calculateMeanFactorToState(const T& sum, const T& factor, const T& mean, const T&  m_state_to_factor){
        const T factorInverse = 1.0/factor;
        return factorInverse * (mean - sum + factor*m_state_to_factor);
    }

    template<typename T>
    static inline const T& calculateVarianceFactorToState(const T& sum, const T& factor, const T& variance, const T& v_state_to_factor){
        const T factorSquared = factor * factor;
        const T factorSquaredInverse = 1.0/factorSquared;
        return factorSquaredInverse * (variance - sum + factorSquared*v_state_to_factor);
    }
};

#endif
