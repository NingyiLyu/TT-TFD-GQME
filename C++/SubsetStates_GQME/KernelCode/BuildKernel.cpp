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
#include "FunctionTemplates_Kernel.h"

// Purpose: To calculate the memory kernel
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created by Ellen Mulvihill on 9/27/18
//
// Most recent modification on 8/29/22

// builds the memory kernel from the non-Condon Volterra equation
void BuildKernel(string& idStr_K, string& idStr_PFI, vector<string>& states,
                 int& numStates){
    int i,j,k,l,m,n; // integers for loops
    
    // liouvillian matrix
    Complex_Matrix LN0(DOF_E_SQ, vector<Complex> (DOF_E_SQ, 0.0));
    MakeExpvN0Liouville(LN0); // in ../../HamiltonianFunctions.cpp

    // calculating the kernel
    cout << "   >>> Building Single-State projection-free input " << endl;
        
    Complex_3D_Matrix linearTerm(LEN_TRAJ_K, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
        
    // The projection-free inputs F and Fdot
    Complex_3D_Matrix F(LEN_TRAJ_K, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
    Complex_3D_Matrix Fdot(LEN_TRAJ_K, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
            
    // reads in the projection-free inputs F and Fdot, function is
    // located in ReadInPFI.cpp
    ReadInPFI(idStr_PFI, states, numStates, "F", F);
    ReadInPFI(idStr_PFI, states, numStates, "Fdot", Fdot);
            
    // creates the linear term of the volterra equation, e.g. the eqn is
    // K(t) = i * Fdot(t) - 1/hbar * F * <L>_n^0 + i * int_0^t F(t - t') K(t')
    // so this calculates i * Fdot(t) - 1/hbar * F * <L>_n^0
    for (i = 0; i < numStates; i++){
        for (j = 0; j < numStates; j++){
            for (l = 0; l < LEN_TRAJ_K; l++){
                linearTerm[l][i][j] = I * Fdot[l][i][j];
                for (k = 0; k < numStates; k++){
                    int index_uv = DOF_E * (states[k].at(0) - '0') + (states[k].at(1) - '0');
                    int index_lm = DOF_E * (states[j].at(0) - '0') + (states[j].at(1) - '0');
                    linearTerm[l][i][j] -= 1./HBAR * F[l][i][k] * LN0[index_uv][index_lm];
                }
            }
        }
    }
    // calculates the volterra equation for the memory kernel
    // K(t) = i * Fdot(t) + i * int_0^t F(t - t') K(t')
    // with the RunVolterra function in IntegralVolterraFunctions.cpp
    cout << "   >>> Starting K Volterra iterative method" << endl;
    RunVolterra(idStr_K, states, numStates, linearTerm, F);
    
    return;
}

