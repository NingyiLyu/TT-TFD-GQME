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
#include "constants.h"
#include "FunctionTemplates_HamNucMode.h"

// Purpose: Nuclear mode functions
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22

// ----------------------------------
// ----- Nuclear Mode Functions -----
// ----------------------------------

// Checks which nuclear mode function to call based on system
void CallInNucModes(int& flag, vector<double>& omega, vector<double>& shifts,
                    vector<double>& popCoeff, vector<double>& cohCoeff){
    if (SYSTEM == "Spin-Boson"){
        flag = ReadInOhmicParams(omega, shifts, popCoeff);
    }
    else {
        cout << "\t ERROR: System not one of the options in CallInNucModes() ";
        cout << "in NuclearModeFunctions.cpp" << endl;
        flag = -1;
    }
    
    return;
}

int ReadInOhmicParams(vector<double>& omega, vector<double>& shifts,
                      vector<double>& popCoeff){
    // for 2-level system
    // read in force field parameters
    string params("");
    stringstream ss;
    ss << fixed << setprecision(1) << "xi" << XI << "wc";
    ss << fixed << setprecision(1) << OMEGA_C;
    ss << "_wmax" << OMEGA_MAX;
    ss << "_dofn" << DOF_N;
    //cout << ss.str() << endl;
    params += ss.str();
    ss.str("");
    ss.clear();
    
    int i,j;
    ifstream infile;
    // omega input
    infile.open("../BathInputs/Ohmic/specdens_parameters/omega_nm_" + params + ".txt");
    if (!infile.is_open()) {
        cout << "\t ERROR: input file omega cannot open"<< endl;
        return -1;
    }
    if (omega.size() != DOF_N) omega.resize(DOF_N,0);
    for (i = 0; i < DOF_N; i++) {
        infile >> omega[i];
    }
    infile.close();
    infile.clear();
    
    // c input
    infile.open("../BathInputs/Ohmic/specdens_parameters/C_nm_" + params + ".txt");
    if (!infile.is_open()) {
        cout << "\t ERROR: input file popCoeff cannot open"<< endl;
        return -1;
    }
    for (i = 0; i < DOF_N; i++) {
        infile >> popCoeff[i];
    }
    infile.close();
    infile.clear();
    
    cout << "   >>> Input force field parameters successfully." << endl;
    
    return 0;
}
