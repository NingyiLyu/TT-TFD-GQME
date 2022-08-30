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

// Purpose: Functions used within the M-GQME code that directly rely on the
//          Hamiltonian and therefore must be changed for every change in the
//          Hamiltonian
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22

// ------------------------------
// ----- All Code Functions -----
// ------------------------------

// Prints information to shell about the system
void PrintSystemInfo(){
    double dt = DT;
    
    cout << "\t GQME type         = " << GQME_TYPE << endl;
    if (GQME_TYPE == "SubsetStates"){
        cout << "\t Subset            = ";
        for (int i = 0; i < STATES.size(); i++){
            cout << STATES[i] << ", ";
        }
        cout << endl;
    }
    else if (GQME_TYPE == "SingleState"){
        cout << "\t State             = " << STATES[0] << endl;
    }
    cout << "\t System            = " << SYSTEM << endl;
    cout << "\t LEN_TRAJ_PFI      = " << LEN_TRAJ_PFI << endl;
    cout << "\t LEN_TRAJ_K        = " << LEN_TRAJ_K << endl;
    cout << "\t LEN_TRAJ_I        = " << LEN_TRAJ_I << endl;
    cout << "\t DOF_N             = " << DOF_N << endl;
    cout << "\t DOF_E             = " << DOF_E << endl;
    cout << "\t DT                = " << dt << endl;
    cout << "\t t_final           = " << dt * LEN_TRAJ_PFI << endl;
    cout << "\t epsilon           = " << EPSILON << endl;
    cout << "\t Gamma_DA          = " << GAMMA_DA << endl;
    cout << "\t beta              = " << BETA << endl;
    cout << "\t Bath Type         = " << BATH_TYPE << endl;
    cout << "\t Bath Sites        = " << BATH_SITES << endl;
    cout << "\t xi                = " << XI << endl;
    cout << "\t omega_c           = " << OMEGA_C << endl;
    cout << "\t omega_max         = " << OMEGA_MAX << endl;
    
    return;
}

// Checks if system matches other parameters
int CheckSystem(){
    if (SYSTEM == "Spin-Boson"){
        if (BATH_SITES == "All" && BATH_TYPE == "Ohmic"){
            return 0;
        }
        else {
            cout << "Spin-Boson does not have matching BATH_SITES or BATH_TYPE";
            cout << " parameter. Expected All and Ohmic." << endl;
            return -1;
        }
    }
    else {
        cout << "System not given in CheckSystem() in";
        cout << "HamiltonianFunctions.cpp." << endl;
        return -1;
    }
    
    return 0;
}

// Creates the ID string
int CreateIDString(stringstream& ss, string type, string& idStr){
    double dt = DT;
    
    ss << "b" << BETA;
    ss << "G" << GAMMA_DA << "_e" << EPSILON;
    ss << "_t" << setprecision(8) << dt;
    ss << "_xi" << XI << "wc" << OMEGA_C << "_wmax" << OMEGA_MAX;
    ss << "_dofn" << DOF_N;
    if (type == "PFI"){
        ss << "_tf" << setprecision(3) << (LEN_TRAJ_PFI - 1) * dt;
    }
    else if (type == "K"){
        ss << "_tf" << setprecision(3) << (LEN_TRAJ_K - 1) * dt;
    }
    else if (type == "I"){
        ss << "_tf" << setprecision(3) << (LEN_TRAJ_I - 1) * dt;
    }
    else {
        cout << "\tERROR: type not PFI, K, or I in CreateIDString. Exit.";
        cout << endl;
        return 1;
    }
    
    idStr += ss.str();
    ss.str("");
    ss.clear();
    
    return 0;
}

// Creates the ID string
int CreateIDString_PFI(string& idStr){
    double dt = DT;
    
    stringstream ss;
    ss.str("");
    ss << "b" << BETA;
    ss << "G" << GAMMA_DA << "_e" << EPSILON;
    ss << "_t" << setprecision(8) << dt;
    ss << "_xi" << XI << "wc" << OMEGA_C << "_wmax" << OMEGA_MAX;
    ss << "_dofn" << DOF_N;
    ss << "_tf" << setprecision(3) << (LEN_TRAJ_PFI - 1) * dt;
    idStr += ss.str();
    ss.str("");
    ss.clear();
    
    return 0;
}

// --------------------------------
// ----- KernelCode Functions -----
// --------------------------------
double KroneckerDelta(int a, int b){
    if (a == b){
        return 1.;
    }
    else {
        return 0.;
    }
}

// -------------------------------------------
// ----- KernelCode & GQMECode Functions -----
// -------------------------------------------

// Creates the <L>_N^0 matrix
void MakeExpvN0Liouville(Complex_Matrix& LN0){
    
    LN0[0][1] = LN0[1][0] = LN0[2][3] = LN0[3][2] = -GAMMA_DA;
    LN0[0][2] = LN0[2][0] = LN0[1][3] = LN0[3][1] = GAMMA_DA;
    LN0[1][1] = 2. * EPSILON;
    LN0[2][2] = -2. * EPSILON;

    return;
}

// ------------------------------
// ----- GQMECode Functions -----
// ------------------------------

// Prints to shell the system info unique to GQME part
void PrintGQMESystemInfo(){
    
    cout << "\t t_mem             = " << MEM_TIME << endl;
    cout << "\t final t           = " << FINAL_TIME << endl;
    cout << "\t CONV_LIMIT        = " << CONV_LIMIT << endl;
    cout << "\t CONV_ALG_BIG_STEP   = " << CONV_ALG_BIG_STEP << endl;
    cout << "\t CONV_ALG_SMALL_STEP = " << CONV_ALG_SMALL_STEP << endl;
    
    return;
}

// Creates the memory time and final time parts of the print string
void CreatePrintString(double& memTime, string& printStr){
    stringstream ss;
    
    ss << "_mt" << memTime << "_finalt" << FINAL_TIME;
    printStr += ss.str();
    ss.str("");
    ss.clear();
    
    return;
}
