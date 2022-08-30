#ifndef _FUNCTION_TEMPLATES_GQME_H_
#define _FUNCTION_TEMPLATES_GQME_H_

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

// Purpose: To define the function templates for the GQME program
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification 8/29/22 (Ellen)

// -------------------------------------------
// ----- Runge-Kutta 4th Order Functions -----
// -------------------------------------------
void RungeKutta_4O(string& idStr_K, string& idStr_I, vector<string>& states,
                   int& numStates);
void PropagateRK4(vector<string>& states, int& numStates, bool& initInSubset,
                  double& time, Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                  Complex_Matrix& iTerm, vector<Complex>& sigma_hold,
                  Complex_Matrix& sigma, double& memTime);
void Calculatef(vector<string>& states, int& numStates, bool& initInSubset,
                double& time, Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                Complex_Matrix& iTerm, Complex_Matrix& sigma,
                vector<Complex>& k, vector<Complex>& f, double& memTime);
void ConvergenceAlgorithm(vector<string>& states, int& numStates,
                          bool& initInSubset, int& initialIndex, double& time,
                          Complex_Matrix& LN0, Complex_3D_Matrix& kernel,
                          Complex_Matrix& iTerm, Complex_Matrix& sigma,
                          double& memTime);
void MemtimeCheckAlgorithm(string& idStr_K, vector<string>& states,
                           int& numStates, bool& initInSubset,
                           int& initialIndex, double& time, Complex_Matrix& LN0,
                           Complex_3D_Matrix& kernel, Complex_Matrix& iTerm,
                           Complex_Matrix& sigma, double& memTime);

// --------------------------------------
// ----- Input and Output Functions -----
// --------------------------------------
void ReadInKernel(string& idStr_K, vector<string>& states, int& numStates,
                  Complex_3D_Matrix& kernel);
void ReadInInhomogeneousTerm(string& idStr_I, vector<string>& states,
                             int& numStates, Complex_Matrix& iTerm);
void Print_sigma(string& idStr_K, vector<string>& states, int& numStates,
                 Complex_Matrix& sigma, double& memTime);
void Print_mainSigma(vector<string>& states, int& numStates, string& outputStr,
                     Complex_Matrix& sigma);

#endif
