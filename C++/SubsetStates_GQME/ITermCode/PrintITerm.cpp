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

// Purpose: To print the elements of the memory kernel
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created by Ellen Mulvihill on 9/26/18
//
// Most recent modification on 8/29/22

// prints the real and imaginary parts of the elements of the memory kernel
void PrintITerm(string& idStr_I, vector<string>& states, int& numStates,
                Complex_Matrix& iTerm){
    int i,j,k,l; // integers for loops
    string outputStr("");
    ofstream outfile; // variable for output file
    
    for (j = 0; j < numStates; j++){
        outputStr = "I_" + states[j] + "_startingIn_" + INITIAL_STATE + "_";

        if (GQME_TYPE == "SubsetStates"){
            outputStr += "subset_";
            for (l = 0; l < numStates; l++){
                outputStr += states[l] + "_";
            }
        }
        outputStr += BATH_TYPE + "_" + DYN_STRING + "_";
            
        // opens the file to print the real part of the kernel
        outfile.open((K_FOLDER + outputStr + idStr_I + ".dat").c_str());
        // prints error message if file does not open
        if (!outfile.is_open()) {
            cout << "\t ERROR: output file of I " << states[j] << "_";
            cout << INITIAL_STATE << " cannot open" << endl;
            cout << K_FOLDER + outputStr + idStr_I + ".dat" << endl;
        }
        // loops to print the time, kernel, and error bars
        for (i = 0 ; i < LEN_TRAJ_I; i++) {
            double ti = DT*i;
            outfile << ti << "\t" << iTerm[i][j].real() << "\t";
            outfile << iTerm[i][j].imag() << endl;
        }
        outfile.close();
        outfile.clear();
    }
    
    return;
}

