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
#include "FunctionTemplates_GQME.h"

// Purpose: To calculate full and reduced-dimensionality GQMEs
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void RungeKutta_4O(string& idStr_K, string& idStr_I, vector<string>& states,
                   int& numStates){
    int i,j,k,l;
    int flag(0); // used for checking initialization
    double time;
    int initialIndex = 0;
    clock_t clockStartGQME;
    clockStartGQME = clock();

    // making the liouvillian matrix
    Complex_Matrix LN0(DOF_E_SQ, vector<Complex>(DOF_E_SQ, 0.));
    MakeExpvN0Liouville(LN0);
    
    Complex_Matrix sigma(FINAL_TIMESTEPS, vector<Complex>(numStates, 0.));
    Complex_3D_Matrix kernel(MEM_TIMESTEPS, Complex_Matrix(numStates, vector<Complex>(numStates, 0.)));
    Complex_Matrix iTerm(I_TIMESTEPS, vector<Complex>(numStates, 0.));

    // Testing to see if the initial state is within the subset of states being
    // calculated and setting the initial index. The initial state is assumed to
    // be of the form rho(0) = rho_n(0) * |alpha><alpha|, so the initial state
    // is always within the Full and PopulationsOnly GQMEs
    bool initInSubset = false;
    if (GQME_TYPE == "Full"){
        initInSubset = true;
        initialIndex = DOF_E * (INITIAL_STATE.at(0) - '0') + (INITIAL_STATE.at(1) - '0');
    }
    else if (GQME_TYPE == "PopulationsOnly"){ //assumes you start in a pop
        initInSubset = true;
        initialIndex = 1 * (INITIAL_STATE.at(0) - '0');
    }
    else if (GQME_TYPE == "SubsetStates"){
        for (i = 0; i < numStates; i++){
            if (states[i] == INITIAL_STATE){
                initialIndex = i;
                initInSubset = true;
            }
        }
    }
    else if (GQME_TYPE == "SingleState"){
        if (states[0] == INITIAL_STATE){
            initialIndex = 0;
            initInSubset = true;
        }
    }
    cout << "   >>> Setting initial density matrix to state " << INITIAL_STATE;
    cout << endl;
    if (initInSubset == true) {
        sigma[0][initialIndex] = 1.;
    }
    
    cout << "   >>> Reading in memory kernel elements" << endl;
    ReadInKernel(idStr_K, states, numStates, kernel);
    
    if (initInSubset == false){
        cout << "   >>> Reading in inhomogeneous term elements" << endl;
        ReadInInhomogeneousTerm(idStr_I, states, numStates, iTerm);
    }
    
    double memTime = 0;
    if (CONV_ALG == true && MEMTIME_CHECK == true){
        cout << "\tERROR: Cannot have both CONV_ALG and MEMTIME_CHECK as true.";
        cout << "Exiting." << endl;
        
        return;
    }
    else if (CONV_ALG == true){
        cout << "   >>> Starting convergence algorithm" << endl;
        ConvergenceAlgorithm(states, numStates, initInSubset, initialIndex,
                             time, LN0, kernel, iTerm, sigma, memTime);
        
        cout << "   >>> Printing sigma" << endl;
        Print_sigma(idStr_K, states, numStates, sigma, memTime);
    }
    else if (MEMTIME_CHECK == true){
        cout << "   >>> Starting memory time check" << endl;
        MemtimeCheckAlgorithm(idStr_K, states, numStates, initInSubset,
                              initialIndex, time, LN0, kernel, iTerm, sigma,
                              memTime);
    }
    else {
        vector<Complex> sigma_hold(numStates, 0.);
        if (initInSubset == true){
            sigma_hold[initialIndex] = 1.;
        }
    
        memTime = MEM_TIME;
        double memTime_print = memTime;
        AdjustTimePrinted(memTime_print);
        cout << "   >>> Starting propagation, memTime = " << memTime_print;
        cout << endl;
        for (l = 0; l < (FINAL_TIMESTEPS - 1); l++){
            time = l * DT;
            
            PropagateRK4(states, numStates, initInSubset, time, LN0, kernel,
                         iTerm, sigma_hold, sigma, memTime);
            
            for (j = 0; j < numStates; j++){
                sigma[l + 1][j] = sigma_hold[j];
            }
        }
        
        cout << "   >>> Printing sigma" << endl;
        Print_sigma(idStr_K, states, numStates, sigma, memTime);
    }
    
    if ((((clock() - clockStartGQME)/(double) CLOCKS_PER_SEC)/60./60.)
        < 2){
        printf("   GQME time: %5.8f minutes.\n",
               ((clock() - clockStartGQME)/(double) CLOCKS_PER_SEC)/60.);
    }
    else {
        printf("   GQME time: %5.8f hours.\n",
        ((clock() - clockStartGQME)/(double) CLOCKS_PER_SEC)/60./60.);
    }
    
    return;
}
