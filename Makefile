COMPILER=g++
CXXFLAGS=-O3 -Wall
EXPERIMENTALFLAGS=-ffast-math 
all: main sparse_matrix belief_propagation WLS
	$(COMPILER) $(CXXFLAGS) build/main.o build/SparseMatrix.o build/BeliefPropagation.o  build/WLS.o -o bin/main 

main: main.cpp
	$(COMPILER) $(CXXFLAGS) -c main.cpp -o build/main.o

sparse_matrix: lib/src/SparseMatrix.cpp
	$(COMPILER) $(CXXFLAGS) -c lib/src/SparseMatrix.cpp -o build/SparseMatrix.o

belief_propagation:	lib/src/BeliefPropagation.cpp
	$(COMPILER) $(CXXFLAGS) -c lib/src/BeliefPropagation.cpp -o build/BeliefPropagation.o

WLS: lib/src/WLS.cpp
	$(COMPILER) $(CXXFLAGS) -c lib/src/WLS.cpp -o build/WLS.o

run: bin/main
	./bin/main $(folder)

performance_bench: performance_bench.cpp
	$(COMPILER) $(CXXFLAGS) -c performance_bench.cpp -o build/performance_bench.o

build_performance: performance_bench sparse_matrix belief_propagation 
	$(COMPILER) $(CXXFLAGS) build/performance_bench.o build/SparseMatrix.o build/BeliefPropagation.o  -o bin/performance 

bench:
	./bin/performance
