#ifndef VECTOR_IO_H
#define VECTOR_IO_H

#include <vector>
#include <string>
#include <fstream>

using namespace std;

class VectorIO {
public:
    template<typename T>
    static vector<T> loadFromFile(string fileName){
        ifstream fileInput(fileName);
        if(fileInput.fail()){
            throw domain_error("Failed opening file while loading vector");
        }
        unsigned height, width;
        fileInput >> height >> width;
        vector<T> vec(height);
        
        for(unsigned i = 0; i < height; i++){
            fileInput >> vec[i];
        }
        return vec;
    }     
};
#endif // VECTOR_IO_H