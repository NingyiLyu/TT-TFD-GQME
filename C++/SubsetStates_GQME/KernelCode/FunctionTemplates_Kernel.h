#ifndef _FUNCTION_TEMPLATES_KERNEL_H_
#define _FUNCTION_TEMPLATES_KERNEL_H_

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
#include "../../FunctionTemplates_HamNucMode.h"

// Purpose: To define the function templates for the GQME Kernel program
//
// Created by Ellen Mulvihill on 11/8/21
//
// Most recent modification 11/9/21

// -------------------------------
// ----- Proj-free Functions -----
// -------------------------------
void ReadInPFI(string& idStr_PFI, vector<string>& states, int& numStates,
               string projFreeStr, Complex_3D_Matrix& projFreeInp);
//      reads in a projection-free input; located in ReadInProjFree.cpp

// ----------------------------
// ----- Kernel Functions -----
// ----------------------------
void BuildKernel(string& idStr_K, string& idStr_PFI, vector<string>& states,
                 int& numStates);
//      builds the memory kernel; located in BuildKernel.cpp
void PrintKernel(string method, string& idStr_K, vector<string>& states,
                 int& numStates, Complex_3D_Matrix& kernel);
//      prints the real and imaginary parts of the memory kernel; located in
//      PrintKernel.cpp


// ------------------------------
// ----- Volterra Functions -----
// ------------------------------
void CalculateIntegral(int& numStates, Complex_3D_Matrix& linearTerm,
                       Complex_3D_Matrix& F, Complex_3D_Matrix& prevKernel,
                       Complex_3D_Matrix& kernel);
//      calculates an integral of the form
//      kernel(t_f) = linearTerm(t_f) [i * Fdot(t_f) - 1/hbar * F(t_f) * Ln0]
//                          + i*\int_(t_0)^(t_f) dt' F(t_f - t') * kernel(t')
//      located in IntegralVolterraFunctions.cpp
void RunVolterra(string& idStr_K, vector<string>& states, int& numStates,
                 Complex_3D_Matrix& linearTerm, Complex_3D_Matrix& F);
//      calculates a Volterra equation of the form
//      kernel(t_f) = linearTerm(t_f) [i * Fdot(t_f) - 1/hbar * F(t_f) * Ln0]
//                          + i*\int_(t_0)^(t_f) dt' F(t_f - t') * kernel(t')
//      iteratively with a set convergence parameter and maximum number of
//      iterations, given in ../constants.h; located in
//      IntegralVolterraFunctions.cpp


// ------------------------------------
// ----- Newton-Raphson Functions -----
// ------------------------------------
void CalculateKernelViaNR(string& idStr_K, vector<string>& states,
                          int& numStates, Complex_3D_Matrix& linearTerm,
                          Complex_3D_Matrix& F);
void CalculateNRMethod(int& numStates, Complex_3D_Matrix& S,
                       Complex_3D_Matrix& initialGuess, string termString,
                       Complex_3D_Matrix& term);
void CalculateConvergence(int& numStates, int& numIter, bool& convFlag,
                          clock_t& clock_start, Complex_3D_Matrix& prevTerm,
                          Complex_3D_Matrix& term);


#endif
