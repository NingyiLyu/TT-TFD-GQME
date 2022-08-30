CXX      := g++
CXXFLAGS := -std=c++11
number = 1
system = Spin-Boson
gqmeType = Single_00
objects = main_Kernel.o BuildKernel.o PrintKernel.o IntegralVolterraFunctions.o NewtonRaphsonFunctions.o ReadInPFI.o

$(system)_$(gqmeType)_main_Kernel$(number).exe: $(objects) NuclearModeFunctions.o HamiltonianFunctions.o
	g++ $(objects) NuclearModeFunctions.o HamiltonianFunctions.o -o $(system)_$(gqmeType)_main_Kernel$(number).exe -std=c++11

$(objects): FunctionTemplates_Kernel.h ../../constants.h ../../gen_constants.h

NuclearModeFunctions.o: ../../NuclearModeFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../NuclearModeFunctions.cpp -o NuclearModeFunctions.o -std=c++11

HamiltonianFunctions.o: ../../HamiltonianFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../HamiltonianFunctions.cpp -o HamiltonianFunctions.o -std=c++11

.PHONY : clean
clean:
	rm -f *.o *.exe
