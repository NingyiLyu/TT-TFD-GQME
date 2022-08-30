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

// Purpose: To read in the projection-free inputs used to build the
//          inhomogeneous term
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created on 7/26/21
//
// Most recent modification on 8/29/22

// reads in a projection-free input F
void ReadInF(string& idStr_PFI, vector<string>& states, int& numStates,
             Complex_3D_Matrix& F){
    int i,j,l,m; // integers for loops
    string indicesStr(""); // string specifiying indices of proj-free
    string inputStr("");
    vector<double> time(LEN_TRAJ_I,0.0); // time vector
    double extra = 0.; // throw-away variable for if there are more columns in
    // the SDBCF input files than will be evaluated in the KernelCode
    ifstream infile; // variable for input files
    
    for (i = 0; i < numStates; i++){
        for (j = 0; j < numStates; j++){
            // vectors for the real and imaginary parts of F
            vector<double> projFreeRe(LEN_TRAJ_I, 0.0);
            vector<double> projFreeIm(LEN_TRAJ_I, 0.0);
            
            inputStr = "F_" + states[i] + states[j] + "_" + BATH_TYPE + "_";
            inputStr += DYN_STRING + "_";
            
            // opens the file with the projection-free input
            infile.open(PFI_FOLDER + inputStr + idStr_PFI + ".dat");
            // prints error message if file does not open
            if (!infile.is_open()) {
                cout << "\t ERROR: re input file of " << DYN_STRING << " ";
                cout << "F " << states[i] << states[j] << " cannot open";
                cout << endl;
                // prints full path + name of file it tried to open, to help
                // find the error
                cout << PFI_FOLDER + inputStr + idStr_PFI + ".dat" << endl;
            }
            // loops pull data from file
            for (l = 0; l < LEN_TRAJ_I; l++){
                infile >> time[l] >> projFreeRe[l] >> projFreeIm[l];
            }
            infile.close();
            infile.clear();
            
            // puts the real and imaginary parts of the projection-free input
            // together
            for (l = 0; l < LEN_TRAJ_I; l++){
                F[l][i][j] = projFreeRe[l] + I * projFreeIm[l];
            }
        }
    }
    
    return;
}

// reads in a projection-free input Z
void ReadInZ(string& idStr_PFI, vector<string>& states, int& numStates,
             Complex_Matrix& Z){
    int i,j,l,m; // integers for loops
    string indicesStr(""); // string specifiying indices of proj-free
    string inputStr("");
    vector<double> time(LEN_TRAJ_I,0.0); // time vector
    double extra = 0.; // throw-away variable for if there are more columns in
    // the SDBCF input files than will be evaluated in the KernelCode
    ifstream infile; // variable for input files
    
    for (i = 0; i < numStates; i++){
        // vectors for the real and imaginary parts of Z
        vector<double> projFreeRe(LEN_TRAJ_I, 0.0);
        vector<double> projFreeIm(LEN_TRAJ_I, 0.0);
    
        inputStr = "F_" + states[i] + INITIAL_STATE + "_" + BATH_TYPE + "_";
        inputStr += DYN_STRING + "_";
            
        // opens the file with the projection-free input
        infile.open(PFI_FOLDER + inputStr + idStr_PFI + ".dat");
        // prints error message if file does not open
        if (!infile.is_open()) {
            cout << "\t ERROR: re input file of " << DYN_STRING << " ";
            cout << "F " << states[i] << INITIAL_STATE << " cannot open for Z";
            cout << endl;
            // prints full path + name of file it tried to open, to help find
            // the error
            cout << PFI_FOLDER + inputStr + idStr_PFI + ".dat" << endl;
        }
        // loops pull data from file
        for (l = 0; l < LEN_TRAJ_I; l++){
            infile >> time[l] >> projFreeRe[l] >> projFreeIm[l];
        }
        infile.close();
        infile.clear();
            
        // puts the real and imaginary parts of the projection-free input
        // together
        for (l = 0; l < LEN_TRAJ_I; l++){
            Z[l][i] = -1. * I * (projFreeRe[l] + I * projFreeIm[l]);
        }
    }
    
    return;
}
