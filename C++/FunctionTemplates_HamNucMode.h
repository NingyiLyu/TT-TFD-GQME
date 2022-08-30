#ifndef _FUNCTION_TEMPLATES_HAMNUCMODE_H_
#define _FUNCTION_TEMPLATES_HAMNUCMODE_H_

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

// Purpose: To define the function templates of the Hamiltonian and Nuclear
// Modes Functions for the M-GQME program
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification 8/29/22

// ------------------------------------------
// ----- HAMILTONIAN FUNCTION TEMPLATES -----
// ------------------------------------------

// ----- All Code Functions -----
// ------------------------------
void PrintSystemInfo();
//      prints information about the system
int CheckSystem();
//      checks if system matches other parameters
int CreateIDString(stringstream& ss, string type, string& idStr);
//      creates the ID string
int CreateIDString_PFI(string& idStr);

// ----- KernelCode Function Templates -----
// -----------------------------------------
double KroneckerDelta(int a, int b);
void MakeExpvN0ZeroLiouville(Complex_Matrix& LzeroN0);
//      creates the <L_zero>_0^N matrix

// ----- KernelCode & GQMECode Function Templates -----
// ----------------------------------------------------
void MakeExpvN0Liouville(Complex_Matrix& LN0);
//      creates the <L>_N^0 matrix

// ----- GQMECode Function Templates -----
// ---------------------------------------
void PrintGQMESystemInfo();
//      prints system info unique to GQME part
void CreatePrintString(double& memTime, string& printStr);
//      creates the memory time and final time parts of the print string

// -------------------------------------------
// ----- NUCLEAR MODE FUNCTION TEMPLATES -----
// -------------------------------------------
void CallInNucModes(int& flag, vector<double>& omega, vector<double>& shifts,
                    vector<double>& popCoeff, vector<double>& cohCoeff);
//       checks which nuclear mode function to call based on system
void CallInNucModesGQME(int& flag, vector<double>& omega,
                        vector<double>& shifts, vector<double>& popCoeff,
                        vector<double>& cohCoeff);
//       checks which nuclear mode function to call based on system in GQMECode
int ReadInOhmicParams(vector<double>& omega, vector<double>& shifts,
                      vector<double>& popCoeff);
//       initializes force field parameters for Ohmic spectral density

#endif
