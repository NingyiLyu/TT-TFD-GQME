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

// Purpose: To propagate the GQME following the fourth-order Runge-Kutta
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void PropagateRK4(vector<string>& states, int& numStates, bool& initInSubset,
                  double& time, Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                  Complex_Matrix& iTerm, vector<Complex>& sigma_hold,
                  Complex_Matrix& sigma, double& memTime){
    int i;
    double t_1, t_2, t_3;
    vector<Complex> f_0(numStates, 0.);
    vector<Complex> f_1(numStates, 0.);
    vector<Complex> f_2(numStates, 0.);
    vector<Complex> f_3(numStates, 0.);
    vector<Complex> k_1(numStates, 0.);
    vector<Complex> k_2(numStates, 0.);
    vector<Complex> k_3(numStates, 0.);
        
    Calculatef(states, numStates, initInSubset, time, LN0, kernel, iTerm,
               sigma, sigma_hold, f_0, memTime);
        
    t_1 = time + DT/2.;
    for (i = 0; i < numStates; i++){
        k_1[i] = sigma_hold[i] + DT * f_0[i]/2.;
    }
    Calculatef(states, numStates, initInSubset, t_1, LN0, kernel, iTerm,
               sigma, k_1, f_1, memTime);
 
    t_2 = time + DT/2.;
    for (i = 0; i < numStates; i++){
        k_2[i] = sigma_hold[i] + DT * f_1[i]/2.;
    }
    
    Calculatef(states, numStates, initInSubset, t_2, LN0, kernel, iTerm,
               sigma, k_2, f_2, memTime);
        
    t_3 = time + DT;
    for (i = 0; i < numStates; i++){
        k_3[i] = sigma_hold[i] + DT * f_2[i];
    }
        
    Calculatef(states, numStates, initInSubset, t_3, LN0, kernel, iTerm,
               sigma, k_3, f_3, memTime);
 
    for (i = 0; i < numStates; i++){
        sigma_hold[i] += DT/6. * (f_0[i] + 2. * f_1[i] + 2. * f_2[i] + f_3[i]);
    }
 
    return;
}

void Calculatef(vector<string>& states, int& numStates, bool& initInSubset,
                double& time, Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                Complex_Matrix& iTerm, Complex_Matrix& sigma,
                vector<Complex>& k, vector<Complex>& f, double& memTime){
    int i,j,l,limit;
    Complex hold = 0.;
    Complex memory = 0.;
    int memTimeSteps = memTime/DT;
    
    for (i = 0; i < numStates; i++){
        f[i] = 0.;
        for (j = 0; j < numStates; j++){
            int index_jk = DOF_E * (states[i].at(0) - '0') + (states[i].at(1) - '0');
            int index_lm = DOF_E * (states[j].at(0) - '0') + (states[j].at(1) - '0');
            f[i] -= I/HBAR * LN0[index_jk][index_lm] * k[j];
        }
        if (initInSubset == false){
            // adding inhomogeneous term if time is less than the limit of the
            // inhomogeneous term
            if (int(time/DT) < I_TIMESTEPS){
                f[i] += iTerm[int(time/DT)][i];
            }
        }
        memory = 0.;
        limit = memTimeSteps;
        if (int(time/DT) < memTimeSteps - 1){
            limit = int(time/DT);
        }
        for (j = 0; j < numStates; j++){
            for (l = 0; l < limit; l++){
                memory -= DT * kernel[l][i][j] * sigma[int(time/DT) - l][j];
            }
        }
        f[i] += memory;
    }

    return;
}
