#ifndef _GEN_CONSTANTS_H_
#define _GEN_CONSTANTS_H_

#include <iostream>
#include <complex>
using namespace std;

// General constants
//
// Created by Ellen Mulvihill 8/29/22
// Most recent modification: 8/29/22 (Ellen)

// -----------------------------------------
// ----- System and Dynamics constants -----
// -----------------------------------------
const string DYN_STRING = "TT-TFD";
const string SYSTEM = "Spin-Boson";
const string GQME_TYPE = "SingleState";
// options: Full, SingleState, SubsetStates, PopulationsOnly
const vector <string> STATES = {"11"};
// The state(s) to be looking at for SingleState or SubsetStates. It isn't
// necessary to set this for Full or PopulationsOnly
const string INITIAL_STATE = "00";
const string TYPE_COUPLING = "Non-Condon";

// -------------------------------------
// ----- Input/output folder names -----
// -------------------------------------
const string OUTPUT_FOLDER = "../../../Output/";
const string SIGMA_FOLDER = OUTPUT_FOLDER + "Sigma_Output/";
const string PFI_FOLDER = OUTPUT_FOLDER + "ProjFree_Output/";

const string GQME_TYPE_FOLDER = GQME_TYPE + "_GQME/";
const string K_FOLDER = OUTPUT_FOLDER + GQME_TYPE_FOLDER + "K_Output/";
const string GQME_FOLDER = OUTPUT_FOLDER + GQME_TYPE_FOLDER + "GQME_Output/";

// -------------------------------------
// ----- Variable type definitions -----
// -------------------------------------
typedef std::complex<double>  Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<Complex> > > Complex_3D_Matrix;
typedef std::vector<vector<vector<double> > > Real_3D_Matrix;
typedef std::vector<vector<vector<vector<Complex> > > > Complex_4D_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4D_Matrix;

// --------------------------------
// ----- Conversion constants -----
// --------------------------------
const double CM_PER_HARTREE = 2.19475e5; // constant for cm^-1/hartree (a.u.)
const double K_PER_HARTREE = 3.15775e5; // constant for K/hartree (a.u.)
const double PS_PER_AU = 2.41888e-5; // constant for ps/a.u.
const double EV_PER_HARTREE = 27.21139; // constant for eV/hartree (a.u.)
const double PS_PER_S = 1.e12; // constant for ps/s

// -----------------------------
// ----- General constants -----
// -----------------------------
const double PI = std::acos(-1.0);
const complex<double> I(0,1);
const double OMEGA_C_UNIT = 1.; //106.14 / CM_PER_HARTREE; // diabatic coupling
//                                for making other params in terms of omega_c;
//                                set to 1 if want a.u. (Dimer/FMO)
const double HBAR = 1. / OMEGA_C_UNIT;
//const double HBAR = 6.582119569e-16 / EV_PER_HARTREE * CM_PER_HARTREE * PS_PER_S;
//                  hbar in terms of cm^{-1} * ps 

#endif
