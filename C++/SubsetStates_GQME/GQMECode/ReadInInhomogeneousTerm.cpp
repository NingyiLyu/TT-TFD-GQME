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

// Purpose: To read in the elements of the inhomogeneous term
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void ReadInInhomogeneousTerm(string& idStr_I, vector<string>& states,
                             int& numStates, Complex_Matrix& iTerm){
    int i,j,l;
    double time_hold;
    string inputStr("");
    ifstream infile;
    
    for (j = 0; j < numStates; j++){
        // vectors for the real and imaginary parts of the inhomogeneous term
        vector<double> iTermRe(I_TIMESTEPS, 0.0);
        vector<double> iTermIm(I_TIMESTEPS, 0.0);
    
        inputStr = "I_" + states[j] + "_startingIn_" + INITIAL_STATE + "_";
        if (GQME_TYPE == "SubsetStates"){
            inputStr += "subset_";
            for (l = 0; l < numStates; l++){
                inputStr += states[l] + "_";
            }
        }
        inputStr += BATH_TYPE + "_" + DYN_STRING + "_" + idStr_I;
        
        infile.open((K_FOLDER + inputStr + ".dat").c_str());
        if (!infile.is_open()) {
            cout << "\t ERROR: re input file of I" << states[j];
            cout << " starting in " <<  INITIAL_STATE << " " << DYN_STRING;
            cout << " cannot open" << endl;
            cout << K_FOLDER + inputStr << ".dat" << endl;
        }
        for (i = 0 ; i < I_TIMESTEPS; i++) {
            infile >> time_hold >> iTermRe[i] >> iTermIm[i];
        }
        infile.close();
        infile.clear();
            
        for (i = 0 ; i < I_TIMESTEPS; i++) {
            iTerm[i][j] = iTermRe[i] + I * iTermIm[i];
        }
    }
    
    return;
}
