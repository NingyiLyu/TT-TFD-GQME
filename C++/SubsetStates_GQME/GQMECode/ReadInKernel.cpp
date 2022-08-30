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

// Purpose: To read in the elements of the memory kernel
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

void ReadInKernel(string& idStr_K, vector<string>& states, int& numStates,
                  Complex_3D_Matrix& kernel){
    int i,j,k,l;
    double time_hold;
    string inputStr("");
    ifstream infile;
    
    for (j = 0; j < numStates; j++){
        for (k = 0; k < numStates; k++){
            // vectors for the real and imaginary parts of the kernel
            vector<double> kernelRe(MEM_TIMESTEPS, 0.0);
            vector<double> kernelIm(MEM_TIMESTEPS, 0.0);
    
            inputStr = "K_" + states[j] + states[k] + "_";
            if (GQME_TYPE == "SubsetStates"){
                inputStr += "subset_";
                for (l = 0; l < numStates; l++){
                    inputStr += states[l] + "_";
                }
            }
            inputStr += BATH_TYPE + "_" + DYN_STRING + "_";

            infile.open((K_FOLDER + inputStr + idStr_K + ".dat").c_str());
            if (!infile.is_open()) {
                cout << "\t ERROR: re input file of K_" << states[j];
                cout << states[k] << " " << DYN_STRING << " cannot open";
                cout << endl;
                cout << K_FOLDER + inputStr + idStr_K << ".dat" << endl;
            }
            for (i = 0 ; i < MEM_TIMESTEPS; i++) {
                infile >> time_hold >> kernelRe[i] >> kernelIm[i];
            }
            infile.close();
            infile.clear();
            
            for (i = 0 ; i < MEM_TIMESTEPS; i++) {
                kernel[i][j][k] = kernelRe[i] + I * kernelIm[i];
            }
        }
    }
    
    return;
}
