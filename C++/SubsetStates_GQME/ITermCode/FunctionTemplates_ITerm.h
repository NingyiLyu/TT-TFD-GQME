#ifndef _FUNCTION_TEMPLATES_ITERM_H_
#define _FUNCTION_TEMPLATES_ITERM_H_

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

// Purpose: To define the function templates for the inhomogeneous term
//
// Created by Ellen Mulvihill on 11/8/21
//
// Most recent modification 11/9/21

// -------------------------------
// ----- Proj-free Functions -----
// -------------------------------
void ReadInF(string& idStr_PFI, vector<string>& states, int& numStates,
             Complex_3D_Matrix& F);
//      reads in F; located in ReadInPFI.cpp
void ReadInZ(string& idStr_PFI, vector<string>& states, int& numStates,
             Complex_Matrix& Z);
//      reads in Z; located in ReadInPFI.cpp

// ----------------------------
// ----- Kernel Functions -----
// ----------------------------
void BuildITerm(string& idStr_I, string& idStr_PFI, vector<string>& states,
                int& numStates);
//      builds the memory kernel; located in BuildKernel.cpp
void PrintITerm(string& idStr_I, vector<string>& states, int& numStates,
                Complex_Matrix& iTerm);
//      prints the real and imaginary parts of the memory kernel; located in
//      PrintKernel.cpp


// ------------------------------
// ----- Volterra Functions -----
// ------------------------------
void CalculateIntegral(int& numStates, Complex_Matrix& Z, Complex_3D_Matrix& F,
                       Complex_Matrix& prevITerm, Complex_Matrix& iTerm);
//      calculates an integral of the form
//      I(t_f) = Z(t_f) + i*\int_(t_0)^(t_f) dt' F(t_f - t') * I(t')
//      located in IntegralVolterraFunctions.cpp
void RunVolterra(string& idStr_I, vector<string>& states, int& numStates,
                 Complex_Matrix& Z, Complex_3D_Matrix& F);
//      calculates a Volterra equation of the form
//      I(t_f) = Z(t_f) + i*\int_(t_0)^(t_f) dt' F(t_f - t') * I(t')
//      iteratively with a set convergence parameter and maximum number of
//      iterations, given in ../constants.h; located in
//      IntegralVolterraFunctions.cpp


#endif
