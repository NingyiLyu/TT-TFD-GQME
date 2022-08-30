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

// Purpose: To calculate the memory time check algorithm of the GQME
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void MemtimeCheckAlgorithm(string& idStr_K, vector<string>& states,
                           int& numStates, bool& initInSubset,
                           int& initialIndex, double& time, Complex_Matrix& LN0,
                           Complex_3D_Matrix& kernel, Complex_Matrix& iTerm,
                           Complex_Matrix& sigma, double& memTime){
    int i,j,k,l;
    clock_t clockStartMemT;
    
    for (memTime = 0; memTime <= MEM_TIME; memTime += MEMTIME_CHECK_STEP){
        clockStartMemT = clock();
        vector<Complex> sigma_hold(numStates, 0.);
        Complex_Matrix sigma_memTime(FINAL_TIMESTEPS, vector<Complex>(numStates, 0.));
       
        if (initInSubset == true){
            sigma_hold[initialIndex] = 1.;
            sigma_memTime[0][initialIndex] = 1.;
        }
        
        double memTime_print = memTime;
        AdjustTimePrinted(memTime_print);
        cout << "   >>> Starting propagation, memTime = " << memTime_print;
        cout << endl;
        for (l = 0; l < (FINAL_TIMESTEPS - 1); l++){
            time = l * DT;
                
            PropagateRK4(states, numStates, initInSubset, time, LN0, kernel,
                         iTerm, sigma_hold, sigma_memTime, memTime);
                
            for (j = 0; j < numStates; j++){
                sigma_memTime[l + 1][j] = sigma_hold[j];
            }
        }
        
        cout << "   >>> Printing sigma" << endl;
        Print_sigma(idStr_K, states, numStates, sigma_memTime, memTime);
        
        if ((((clock() - clockStartMemT)/(double) CLOCKS_PER_SEC)/60./60.)
            < 2){
            printf("   \tMemory time: %5.8f minutes.\n",
                   ((clock() - clockStartMemT)/(double) CLOCKS_PER_SEC)/60.);
        }
        else {
            printf("   \tMemory time: %5.8f hours.\n",
            ((clock() - clockStartMemT)/(double) CLOCKS_PER_SEC)/60./60.);
        }
    }
    
    return;
}
