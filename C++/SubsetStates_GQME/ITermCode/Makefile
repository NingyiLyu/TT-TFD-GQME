CXX      := g++
CXXFLAGS := -std=c++11
number = 1
gqmeType = Single_11
objects = main_ITerm.o BuildITerm.o PrintITerm.o IntegralVolterraFunctions.o ReadInPFI.o

$(gqmeType)_main_ITerm$(number).exe: $(objects) NuclearModeFunctions.o HamiltonianFunctions.o ../../constants.h ../../gen_constants.h
	g++ $(objects) NuclearModeFunctions.o HamiltonianFunctions.o -o $(gqmeType)_main_ITerm$(number).exe -std=c++11

$(objects): FunctionTemplates_ITerm.h ../../constants.h ../../gen_constants.h

NuclearModeFunctions.o: ../../NuclearModeFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../NuclearModeFunctions.cpp -o NuclearModeFunctions.o -std=c++11

HamiltonianFunctions.o: ../../HamiltonianFunctions.cpp ../../constants.h ../../gen_constants.h
	g++ -c ../../HamiltonianFunctions.cpp -o HamiltonianFunctions.o -std=c++11

.PHONY : clean
clean:
	rm -f *.o *.exe
