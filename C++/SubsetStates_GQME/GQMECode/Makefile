CXX      := g++
CXXFLAGS := -std=c++11
number = 1
gqmeType = Single_11
objects = main_GQME.o PrintFunctions.o ReadInKernel.o ReadInInhomogeneousTerm.o RungeKutta_4O.o ConvergenceAlgorithm.o MemtimeCheckAlgorithm.o PropagateRK4.o 

$(gqmeType)_main_GQME$(number).exe: $(objects) NuclearModeFunctions.o HamiltonianFunctions.o
	g++ $(objects) NuclearModeFunctions.o HamiltonianFunctions.o -o $(gqmeType)_main_GQME$(number).exe -std=c++11

$(objects): FunctionTemplates_GQME.h ../../constants.h ../../gen_constants.h

NuclearModeFunctions.o: ../../NuclearModeFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../NuclearModeFunctions.cpp -o NuclearModeFunctions.o -std=c++11

HamiltonianFunctions.o: ../../HamiltonianFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../HamiltonianFunctions.cpp -o HamiltonianFunctions.o -std=c++11

.PHONY : clean
clean:
	rm -f *.o *.exe
