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

// Purpose: To print the elements of the memory kernel
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created by Ellen Mulvihill on 9/26/18
//
// Most recent modification on 8/29/22

// prints the real and imaginary parts of the elements of the memory kernel
void PrintKernel(string method, string& idStr_K, vector<string>& states,
                 int& numStates, Complex_3D_Matrix& kernel){
    int i,j,k,l; // integers for loops
    string outputStr("");
    ofstream outfile; // variable for output file
    
    for (j = 0; j < numStates; j++){
        for (k = 0; k < numStates; k++){
            outputStr = "K_" + method + states[j] + states[k] + "_";
            if (GQME_TYPE == "SubsetStates"){
                outputStr += "subset_";
                for (l = 0; l < numStates; l++){
                    outputStr += states[l] + "_";
                }
            }
            outputStr += BATH_TYPE + "_" + DYN_STRING + "_";
            
            // opens the file to print the real part of the kernel
            outfile.open((K_FOLDER + outputStr + idStr_K + ".dat").c_str());
            // prints error message if file does not open
            if (!outfile.is_open()) {
                cout << "\t ERROR: output file of K " << states[j] << states[k];
                cout << " cannot open" << endl;
                cout << K_FOLDER + outputStr + idStr_K + ".dat" << endl;
            }
            // loops to print the time, kernel, and error bars
            for (i = 0; i < LEN_TRAJ_K; i++) {
                double ti = DT*i;
                outfile << ti << "\t" << kernel[i][j][k].real() << "\t";
                outfile << kernel[i][j][k].imag() << endl;
            }
            outfile.close();
            outfile.clear();
        }
    }
    
    return;
}

