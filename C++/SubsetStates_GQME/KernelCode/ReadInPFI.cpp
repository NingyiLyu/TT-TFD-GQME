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

// Purpose: To read in the projection-free inputs used to build the memory
//          kernel from files created by DynamicsCode and SDBCFCode
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22

// reads in a projection-free input
void ReadInPFI(string& idStr_PFI, vector<string>& states, int& numStates,
               string projFreeStr, Complex_3D_Matrix& projFreeInp){
    int i,j,l,m; // integers for loops
    string inputStr("");
    vector<double> time(LEN_TRAJ_K,0.0); // time vector
    double extra = 0.; // throw-away variable for if there are more columns in
    // the SDBCF input files than will be evaluated in the KernelCode
    ifstream infile; // variable for input files
    
    for (i = 0; i < numStates; i++){
        for (j = 0; j < numStates; j++){
            // vectors for the real and imaginary parts of the PFIs
            vector<double> projFreeRe(LEN_TRAJ_K, 0.0);
            vector<double> projFreeIm(LEN_TRAJ_K, 0.0);
            
            //inputStr = projFreeStr + "_blockBLOckNUM_" + states[i] + states[j] + "_";
            inputStr = projFreeStr + "_" + states[i] + states[j] + "_";
            inputStr += BATH_TYPE + "_" + DYN_STRING + "_";
            
            // opens the file with the projection-free input
            infile.open(PFI_FOLDER + inputStr + idStr_PFI + ".dat");
            // prints error message if file does not open
            if (!infile.is_open()) {
                cout << "\t ERROR: re input file of " << DYN_STRING << " ";
                cout << projFreeStr << " " << states[i] << states[j];
                cout << " cannot open" << endl;
                // prints full path + name of file it tried to open, to help find the
                // error
                cout << PFI_FOLDER + inputStr + idStr_PFI + ".dat" << endl;
            }
            // loops pull data from file
            for (l = 0; l < LEN_TRAJ_K; l++){
                infile >> time[l] >> projFreeRe[l] >> projFreeIm[l];
            }
            infile.close();
            infile.clear();
            
            // puts the real and imaginary parts of the projection-free input
            // together
            for (l = 0; l < LEN_TRAJ_K; l++){
                projFreeInp[l][i][j] = projFreeRe[l] + I * projFreeIm[l];
            }
        }
    }
    
    return;
}
