#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "SparseMatrix.hpp"
#include <algorithm>
#include <map>
class Preprocessing {
        static constexpr double virtualMean = 0;
        static constexpr double virtualVariance = 1e3;
public:
    static constexpr double lower_limit = 1e-60, upper_limit = 1e60;

    template<typename T>
    static SparseMatrix<T> reorderNodes(SparseMatrix<T> &factors, vector<T> &means, vector<T> &variances){
        vector<vector<Node<T>>> data = factors.getData();
        multimap<unsigned, int> index_map;
        for(unsigned i = 0; i < data.size(); i++){
            if(data[i].size() == 1){
                index_map.emplace(data[i][0].index, i);
            }
        }
        unsigned unique_single = 0;
        for(auto i = index_map.begin(); i != index_map.end(); ){
            unique_single++;
            i = index_map.upper_bound(i->first);
        }
        unsigned new_size = factors.getHeight() - index_map.size() + unique_single;
        vector<vector<Node<T>>> new_data(new_size);
        vector<T> new_means(new_size), new_varinaces(new_size);

        unsigned counter = 0;
        for(auto i = index_map.begin(); i != index_map.end(); ){
            
            auto range = index_map.equal_range(i->first);
            unsigned index;
            T mean = T(), variance = T();
            for(auto j = range.first; j != range.second; j++){
                index = j->second;
                T new_mean = means[index]/data[index][0].value;
                T new_variance = variances[index]/(data[index][0].value*data[index][0].value);
                if(j == range.first) {
                    mean = new_mean/new_variance;
                    variance = 1.0/new_variance;
                } else {
                    mean += new_mean/new_variance;
                    variance += 1.0/new_variance;
                }
            }
            new_data[counter].emplace_back(Node<T>{1, i->first});
            new_means[counter] = mean/variance;
            new_varinaces[counter] = 1/variance;
            
            counter++;
            i = range.second;
        }
        
        for(unsigned i = 0; i < factors.getHeight(); i++){
            if(data[i].size() != 1) {
                swap(new_data[counter], data[i]);
                new_means[counter] = means[i];
                new_varinaces[counter] = variances[i];
                counter++;
            }
        }
        swap(means, new_means);
        swap(variances, new_varinaces);
        return SparseMatrix<T>(new_data.size(), factors.getWidth(), new_data);
    }

    template<typename T>
    static SparseMatrix<T> addVirtualNodes(SparseMatrix<T> &factors, vector<T>& means, vector<T>& variances){
      
        vector<vector<Node<T>>> data = factors.getData();
        unsigned single_connected_count = 0;
        while(factors.getRowSize(single_connected_count) == 1){
            single_connected_count++;
        }
        unsigned new_height = factors.getHeight() - single_connected_count + factors.getWidth();
        vector<vector<Node<T>>> new_data(new_height);
        vector<T> new_means(new_height), new_variances(new_height);
        unsigned counter = 0;

        for(unsigned i = 0; i < single_connected_count; i++){
            while(counter != data[i][0].index){
                new_data[counter].push_back(Node<T>({1, counter}));
                new_means[counter] = virtualMean;
                new_variances[counter] = virtualVariance;
                counter++;
            }
            new_data[counter].push_back(Node<T>({1, counter}));
            new_means[counter] = means[i];
            new_variances[counter] = variances[i];
            counter++;
        }
        for(unsigned i = single_connected_count; i < data.size(); i++){
            swap(new_data[counter], data[i]);
            new_means[counter] = means[i];
            new_variances[counter] = variances[i];
            counter++; 
        }
        swap(new_means, means);
        swap(new_variances, variances);
        return SparseMatrix<T>(new_data.size(), factors.getWidth(), new_data);
    }

    template<typename T>
    static void normalizeVariances(vector<T> &variances){
        for(unsigned i = 0; i < variances.size(); i++){
            if(variances[i] < lower_limit)
                variances[i] = lower_limit;
            if(variances[i] > upper_limit)
                variances[i] = upper_limit;
        }
    }
    template<typename T>
    static void removeDirectNodes(SparseMatrix<T> &factors, vector<T> &means, vector<T> &variances, vector<T> &direct_means, vector<T> &direct_variances){
        vector<bool> is_connected(factors.getWidth(), false);
        direct_means = vector<T>(factors.getWidth(), 0);
        direct_variances = vector<T>(factors.getWidth(), 0);
        vector<vector<Node<T>>> data = factors.getData();
        unsigned direct_connected = 0;
        for(unsigned i = 0; i < data.size(); i++){
            if(data[i].size() == 1){
                unsigned index = data[i][0].index;
                const T coeff = data[i][0].value;

                direct_means[index] += coeff*means[i]/variances[i];
                direct_variances[index] += coeff*coeff/variances[i];
                is_connected[index] = true;
                direct_connected++;
            }
        }
        for(unsigned i = 0; i < direct_means.size(); i++){
            if(is_connected[i] == false){
                direct_means[i] = virtualMean/virtualVariance;
                direct_variances[i] = 1/virtualVariance;
            }
        }

        unsigned counter = 0, new_height = factors.getHeight() - direct_connected;
        vector<vector<Node<T>>> new_data(new_height);
        vector<T> new_means(new_height), new_variances(new_height);
        for(unsigned i = 0; i < data.size(); i++){
            if(data[i].size() != 1){
                swap(new_data[counter], data[i]);
                new_means[counter] = means[i];
                new_variances[counter] = variances[i];
                counter++;
            }
        }
        swap(means, new_means);
        swap(variances, new_variances);
        factors = SparseMatrix<T>(new_height, factors.getWidth(), new_data);
    }
};

#endif // PREPROCESSING_H