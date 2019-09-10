COMPILER=g++
FLAGS=-O3 -Wall

all: main sparse_matrix belief_propagation  
	$(COMPILER) $(FLAGS) build/main.o build/SparseMatrix.o build/BeliefPropagation.o  -o bin/main 

main: main.cpp
	$(COMPILER) $(FLAGS) -c main.cpp -o build/main.o

sparse_matrix: lib/src/SparseMatrix.cpp
	$(COMPILER) $(FLAGS) -c lib/src/SparseMatrix.cpp -o build/SparseMatrix.o

belief_propagation:	lib/src/BeliefPropagation.cpp
	$(COMPILER) $(FLAGS) -c lib/src/BeliefPropagation.cpp -o build/BeliefPropagation.o

run: bin/main
	./bin/main $(folder)