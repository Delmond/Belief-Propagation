#ifndef ANALYTICS_H
#define ANALYTICS_H
#include "SparseMatrix.hpp"

class Analytics {
    public:
    template<typename T>
    static double computeWrss(vector<T> marginals, const SparseMatrix<T> &H, vector<T> means, vector<T> variances){
        double wrss = 0.0;

        #pragma omp parallel for reduction(+:wrss)
        for(unsigned i = 0; i < H.getHeight(); i++){
            double sum = 0.0;
            for(unsigned j = 0; j < H.getRowSize(i); j++){
                sum += H(i,j) * marginals[H.getIndex(i, j)];
            }
            wrss += (sum - means[i])*(sum - means[i])/variances[i];
        }
        return wrss;
    }
};

#endif // ANALYTICS_H