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

// Purpose: To print the propagated reduced density matrix and the smaller
//          blocks of trajectories for error bars for full GQME dynamics
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification 8/29/22 (Ellen)

void Print_sigma(string& idStr_K, vector<string>& states, int& numStates,
                 Complex_Matrix& sigma, double& memTime){
    int i,j,k, errNum, stepSizeIncrement;
    string outputStr;
    string printStr("");
    
    printStr = DYN_STRING + "_" + idStr_K;
    CreatePrintString(memTime, printStr);
    
    outputStr = BATH_TYPE + "_RK4_" + printStr + ".dat";
    
    Print_mainSigma(states, numStates, outputStr, stepSizeIncrement, sigma);
    
    return;
}

void Print_mainSigma(vector<string>& states, int& numStates, string& outputStr,
                     Complex_Matrix& sigma){
    int i,j,l;
    string conv("");
    string outputStr2("");
    ofstream outfile;
    
    if (CONV_ALG == true){
        conv = "conv_";
    }
    
    for (l = 0; l < numStates; l++){
        outputStr2 = "Sigma_" + states[l] + "_startingIn_" + INITIAL_STATE + "_";
        
        if (GQME_TYPE == "SubsetStates"){
            outputStr2 += "subset_";
            for (j = 0; j < numStates; j++){
                outputStr2 += states[j] + "_";
            }
        }
        outfile.open((GQME_FOLDER + conv + outputStr2 + outputStr).c_str());
        if (!outfile.is_open()) {
            cout << "\t ERROR: output file of sigma " << states[l] << "_";
            cout << INITIAL_STATE << " cannot open" << endl;
            cout << GQME_FOLDER + "Sigma_" + states[l] + "_startingIn_";
            cout << INITIAL_STATE + "_" + outputStr << endl;
        }
        for (i = 0 ; i < FINAL_TIMESTEPS; i++) {
            double ti = DT * i;
            AdjustTimePrinted(ti);
            outfile << ti << "\t" << sigma[i][l].real() << "\t";
            outfile << sigma[i][l].imag() << endl;
        }
        outfile.close();
        outfile.clear();
    }
    
    return;
}
