#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>
#include <complex>
using namespace std;
#include "../../constants.h"
#include "FunctionTemplates_ITerm.h"

// Purpose: To calculate the memory kernel of the various methods
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created by Ellen Mulvihill on 9/27/18
//
// Most recent modification on 8/29/22

// builds the inhomogeneous from the non-Condon Volterra equation
void BuildITerm(string& idStr_I, string& idStr_PFI, vector<string>& states,
                int& numStates){
    int i,j,k,l,m,n; // integers for loops

    // calculating the kernel
    cout << "   >>> Reading in scalar projection-free input ";
    cout << endl;
        
    // calculates the Condon memory kernel
    Complex_3D_Matrix F(LEN_TRAJ_I, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
    Complex_Matrix Z(LEN_TRAJ_I, vector<Complex>(numStates, 0.0));
            
    // reads in the projection-free inputs F and Fdot, function is
    // located in ReadInProjFree.cpp
    ReadInF(idStr_PFI, states, numStates, F);
    ReadInZ(idStr_PFI, states, numStates, Z);
    
    // calculates the volterra equation for the memory kernel
    // I(t) = Z(t) + i * F(t) * sigma_{alpha alpha}(0) --> usually zero
    //          + i * int_0^t F(t - t') I(t')
    // with the RunVolterra function in IntegralVolterraFunctions.cpp
    cout << "   >>> Starting scalar I Volterra " << endl;
    RunVolterra(idStr_I, states, numStates, Z, F);
    
    return;
}

