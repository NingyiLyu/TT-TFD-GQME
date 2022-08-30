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

// Purpose: To calculate inhomogeneous term from projection-free inputs
//
// To compile: make
// To execute: ./main_ITerm{number}.exe
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created on 7/20/21
//
// Most recent modification on 8/29/22

int main(int argc, char* argv[]) {
    int i,j,k; // integers for loops
    int flag(0); // used for checking initialization
    stringstream ss; //generic stringstream used for creating ID strings
    string idStr_I(""); //ID string for I with parameters to differentiate files
    string idStr_PFI(""); //ID string for PFI with parameters to differentiate files
    vector<string> states;
    int numStates;
    
    //starting timing
    clock_t clockStart;
    clockStart = clock();
    
    // Creates the states vector and defines the number of states
    if (GQME_TYPE == "Full"){
        numStates = DOF_E_SQ;
        for (i = 0; i < DOF_E; i++){
            for (j = 0; j < DOF_E; j++){
                string statesStr = to_string(i) + to_string(j);
                states.push_back(statesStr);
            }
        }
    }
    else if (GQME_TYPE == "PopulationsOnly"){
        numStates = DOF_E;
        for (i = 0; i < DOF_E; i++){
            string statesStr = to_string(i) + to_string(i);
            states.push_back(statesStr);
        }
    }
    else if (GQME_TYPE == "SubsetStates"){
        numStates = STATES.size();
        for (i = 0; i < STATES.size(); i++){
            states.push_back(STATES[i]);
        }
    }
    else if (GQME_TYPE == "SingleState"){
        numStates = STATES.size();
        if (numStates != 1) {
            cout << "\t ERROR: More than one state in STATES with GQME_TYPE = ";
            cout << "SingleState. Exiting." << endl;
            return -1;
        }
        for (i = 0; i < STATES.size(); i++){
            states.push_back(STATES[i]);
        }
    }
    else {
        cout << "\t ERROR: GQME_TYPE not Full, PopulationsOnly, SubsetStates, ";
        cout << "or SingleState. Exiting." << endl;
        return -1;
    }
    
    // Output to shell to confirm system information
    cout << "   Calculating inhomogeneous term of the " << GQME_TYPE << " GQME";
    cout << endl;
    cout << "   Initial state is " << INITIAL_STATE << endl;
    if (GQME_TYPE == "SubsetStates"){
        cout << "   Subset of states is ";
        for (i = 0; i < numStates; i++){
            cout << states[i] << ", ";
        }
        cout << endl;
    }
    else if (GQME_TYPE == "SingleState"){
        cout << "   State being calculated is " << states[0] << endl;
    }
    cout << "   Dynamics method is " << DYN_STRING << endl;
    
    // Output to the shell, to confirm the constant variables
    PrintSystemInfo();
    
    // Check that parameters match the chosen system
    flag = CheckSystem();
    // Terminates program if error in matching system to params
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in system and parameter ";
        cout << "matching." << endl;
        return -1;
    }
    
    //condition information for ID string
    flag = CreateIDString(ss, "I", idStr_I);
    // Terminates program if error in creating ID string
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in creating I ID string.";
        cout << endl;
        return -1;
    }
    
    //condition information for SDBCF ID string
    flag = CreateIDString(ss, "PFI", idStr_PFI);
    // Terminates program if error in creating ID string
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in creating PF ID string.";
        cout << endl;
        return -1;
    }

    // Calls the BuildKernel function located in BuildKernel.cpp
    BuildITerm(idStr_I, idStr_PFI, states, numStates);
    
    // Prints the total time of the program in minutes if time < 2 hours or in
    // hours if time > 2 hours
    if (((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60. < 120.){
        printf("   Runtime: %5.8f minutes.\n",
               ((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60.);
    }
    else {
        printf("   Runtime: %5.8f hours.\n",
               ((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60./60.);
    }
    
    return 0;
}
