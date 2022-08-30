#include <iostream>
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

// Purpose: To calculate full and reduced-dimensionality GQMEs
//
// To compile: make
// To execute: ./main_GQME{number}.exe
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22 (Ellen)

int main(int argc, char* argv[]){
    int i,j,k,l;
    int flag(0); // used for checking initialization
    stringstream ss; //generic stringstream used for creating ID strings
    string idStr_K(""); //ID string with parameters to differentiate files
    string idStr_I(""); //ID string with parameters to differentiate files
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
    
    cout << "   Calculating scalar GQME type " << GQME_TYPE << endl;
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
    //cout << "   RK4 propagation type is " << typePropRK4 << endl;
    cout << "   Dynamics method is " << DYN_STRING << endl;
    
    // Output to the command line, to confirm the constant variables
    PrintSystemInfo();
    PrintGQMESystemInfo();
    
    // System check
    flag = CheckSystem();
    // Terminates program if error in matching system to params
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in system and parameter ";
        cout << "matching." << endl;
        return -1;
    }
    
    //condition information for ID string
    flag = CreateIDString(ss, "K", idStr_K);
    // Terminates program if error in creating ID string
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in creating K ID string.";
        cout << endl;
        return -1;
    }
    
    flag = CreateIDString(ss, "I", idStr_I);
    // Terminates program if error in creating ID string
    if (flag != 0) {
        cout << "\t ERROR: Terminate due to error in creating K ID string.";
        cout << endl;
        return -1;
    }
    
    if (MEM_TIMESTEPS > LEN_TRAJ_K){
        cout << "\t ERROR: memory time steps more than length of trajectory ";
        cout << "for kernel. Exiting." << endl;
        
        return -1;
    }
    if (I_TIMESTEPS > LEN_TRAJ_I){
        cout << "\t ERROR: I time steps more than length of trajectory for ";
        cout << "inhomogeneous term. Exiting." << endl;
        
        return -1;
    }
    
    cout << "   >>> Starting RK4 for population-only GQME " << DYN_STRING;
    cout << endl;
    RungeKutta_4O(idStr_K, idStr_I, states, numStates);
    
    if ((((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60./60.) < 2){
        printf("   Runtime: %5.8f minutes.\n",
               ((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60.);
    }
    else {
        printf("   Runtime: %5.8f hours.\n",
        ((clock() - clockStart)/(double) CLOCKS_PER_SEC)/60./60.);
    }
    
    return 0;
}
