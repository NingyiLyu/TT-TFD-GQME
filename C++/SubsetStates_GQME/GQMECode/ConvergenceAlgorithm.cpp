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

// Purpose: To calculate the convergence algorithm of the GQME
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void ConvergenceAlgorithm(vector<string>& states, int& numStates,
                          bool& initInSubset, int& initialIndex, double& time,
                          Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                          Complex_Matrix& iTerm, Complex_Matrix& sigma,
                          double& memTime){
    int i,j,k,l, numLoop, convParam;
    double memTime_print = 0.;
    bool convFirst = false;
    bool convSecond = false;
    bool convThird = false;
    clock_t clockStartInit;
    clockStartInit = clock();
    Complex_Matrix sigma_full(FINAL_TIMESTEPS, vector<Complex>(numStates, 0.));
    Complex_Matrix sigma_oneBack(FINAL_TIMESTEPS, vector<Complex>(numStates, 0.));
    vector<Complex> sigma_hold(numStates, 0.);
    
    if (initInSubset == true){
        sigma_full[0][initialIndex] = 1.;
        sigma_hold[initialIndex] = 1.;
    }
    
    memTime = MEM_TIME;
    memTime_print = memTime;
    AdjustTimePrinted(memTime_print);
    cout << "   >>> Starting propagation, memTime = " << memTime_print << endl;
    
    for (l = 0; l < (FINAL_TIMESTEPS - 1); l++){
        time = l * DT;
            
        PropagateRK4(states, numStates, initInSubset, time, LN0, kernel,
                     iTerm, sigma_hold, sigma_full, memTime);
            
        for (i = 0; i < numStates; i++){
            sigma_full[l + 1][i] = sigma_hold[i];
        }
    }
    printf(" \t Time of initial propagation: %5.8f minutes.\n",
           ((clock() - clockStartInit)/(double) CLOCKS_PER_SEC)/60.);

    memTime = MEM_TIME/2.;
    for (numLoop = 0; numLoop < MAX_LOOPS_CONV; numLoop++){
        clock_t clockStartIter;
        clockStartIter = clock();
        Complex_Matrix sigma_iter(FINAL_TIMESTEPS, vector<Complex>(numStates, 0.));
        
        for (i = 0; i < numStates; i++){
            sigma_hold[i] = 0.;
        }
        
        if (initInSubset == true){
            sigma_hold[initialIndex] = 1.;
            sigma_iter[0][initialIndex] = 1.;
        }
        
        memTime_print = memTime;
        AdjustTimePrinted(memTime_print);
        cout << "   >>> Starting propagation, memTime = " << memTime_print;
        cout << endl;
        for (l = 0; l < (FINAL_TIMESTEPS - 1); l++){
            time = l * DT;
                
            PropagateRK4(states, numStates, initInSubset, time, LN0, kernel,
                         iTerm, sigma_hold, sigma_iter, memTime);
                
            for (i = 0; i < numStates; i++){
                sigma_iter[l + 1][i] = sigma_hold[i];
            }
        }
        
        convParam = 0;
        for (j = 0; j < numStates; j++){
            for (i = 0; i < MEM_TIMESTEPS; i++){
                if (abs(sigma_iter[i][j] - sigma_full[i][j]) <= CONV_LIMIT){
                    convParam += 1;
                }
            }
        }
        if (convFirst == false && convParam != (MEM_TIMESTEPS * numStates)){
            cout << "\t Convergence not found at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            memTime += CONV_ALG_BIG_STEP;
        }
        else if (convFirst == false
                 && convParam == (MEM_TIMESTEPS * numStates)){
            cout << "\t First convergence found at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            convFirst = true;
            for (j = 0; j < numStates; j++){
                for (i = 0; i < FINAL_TIMESTEPS; i++){
                    sigma_oneBack[i][j] = sigma_iter[i][j];
                }
            }
            memTime -= CONV_ALG_BIG_STEP;
        }
        else if (convFirst == true && convParam == (MEM_TIMESTEPS * numStates)
                 && convSecond == false){
            cout << "\t Convergence still true at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            memTime -= CONV_ALG_BIG_STEP;
            for (j = 0; j < numStates; j++){
                for (i = 0; i < FINAL_TIMESTEPS; i++){
                    sigma_oneBack[i][j] = sigma_iter[i][j];
                }
            }
        }
        else if (convFirst == true && convParam != (MEM_TIMESTEPS * numStates)
                 && convSecond == false){
            cout << "\t Convergence now false at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            memTime += CONV_ALG_BIG_STEP - CONV_ALG_SMALL_STEP;
            convSecond = true;
        }
        else if (convFirst == true && convParam == (MEM_TIMESTEPS * numStates)
                 && convSecond == true){
            cout << "\t Convergence still true at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            memTime -= CONV_ALG_SMALL_STEP;
            for (j = 0; j < numStates; j++){
                for (i = 0; i < FINAL_TIMESTEPS; i++){
                    sigma_oneBack[i][j] = sigma_iter[i][j];
                }
            }
        }
        else if (convFirst == true && convParam != (MEM_TIMESTEPS * numStates)
                 && convSecond == true){
            cout << "\t Convergence now false at memTime = " << memTime_print;
            cout << endl;
            printf(" \t Time of iteration: %5.8f minutes.\n",
                   ((clock() - clockStartIter)/(double) CLOCKS_PER_SEC)/60.);
            convThird = true;
            memTime += CONV_ALG_SMALL_STEP;
            for (j = 0; j < numStates; j++){
                for (i = 0; i < FINAL_TIMESTEPS; i++){
                    sigma[i][j] = sigma_oneBack[i][j];
                }
            }
        }
        if (convThird == true){
            return;
        }
        else if (numLoop == MAX_LOOPS_CONV - 1){
            cout << "\t ERROR: Convergence not reached." << endl;
            return;
        }
    }
    
    return;
}
