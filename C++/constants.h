#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <iostream>
#include <complex>
#include "gen_constants.h"
using namespace std;

// Constants for the M-GQME program
//
// Created by Ellen Mulvihill on 8/29/22
//
// Most recent modification on 8/29/22

// ----------------------------
// ----- System constants -----
// ----------------------------
const int LEN_TRAJ_PFI = 10200; // steps in on trajectory for PFIs
const int LEN_TRAJ_K = LEN_TRAJ_PFI; // steps in on trajectory for K
const int LEN_TRAJ_I = LEN_TRAJ_K; // steps in on trajectory for I
const double EPSILON = 1; // half of energy splitting
const double BETA = 5;
// inverse finite temperature (Spin-Boson)
const double GAMMA_DA = 1.; // diabatic coupling (Condon) (Spin-Boson)
const double DT = 0.00150082999505279;
const int DOF_N = 60; // nuclear DOF //96 MIA, 78 BMA
const int DOF_E = 2; // electronic DOF
const int DOF_E_SQ = DOF_E * DOF_E; // electronic DOF squared

// ----------------------------------
// ----- Nuclear mode constants -----
// ----------------------------------
const string BATH_SITES = "All"; // whether there is a bath for each
// electronic state or one for all and if the ones for each electronic state
// are the same or different; should be All for spin-boson
const string BATH_TYPE = "Ohmic"; // type of bath
// For Spin-Boson, BATH_TYPE should be Ohmic
// Ohmic nuclear modes
const double XI = 0.1; // friction coefficient, Ohmic
const double OMEGA_C = 1; // cutoff frequency
const int OMEGA_MAX = 5; // maximum frequency for spec dens

// ---------------------------------
// ----- Kernel Code constants -----
// ---------------------------------
const double INTEGRAL_STEP_K = (LEN_TRAJ_K * DT - DT)/(LEN_TRAJ_K - 1.); 
// size of slice in integral, often written h = (x_b - x_a)/N
// [[Note: N = LEN_TRAJ - 1 because there are LEN_TRAJ - 1 slices, since 
// LEN_TRAJ includes a count for t = 0]]
const double INTEGRAL_STEP_I = (LEN_TRAJ_I * DT - DT)/(LEN_TRAJ_I - 1.);
// size of slice in integral, often written h = (x_b - x_a)/N
// [[Note: N = LEN_TRAJ - 1 because there are LEN_TRAJ - 1 slices, since 
// LEN_TRAJ includes a count for t = 0]]
const int MAX_ITERS = 30; // maximum number of iterations for K
const double CONVERGENCE_PARAM = pow(10., -10); // convergence parameter for
// KernelCode volterra iterations

// -------------------------------
// ----- GQME Code constants -----
// -------------------------------
const bool CONV_ALG = false;
// whether or not to go through the convergence
// algorithm. If true, the conv. alg. with calculate sigma(t) of the mem time
// below and use it as the sigma_z(t) that needs to be converged to, referred
// to as sigma_{z,max}(t). If false,the program will calculate sigma(t) for the
// mem time below
const bool MEMTIME_CHECK = false;
const double MEMTIME_CHECK_START = 2.;//MEM_TIME * 2/3.;
const double MEMTIME_CHECK_STEP = 2.;//(MEM_TIME - MEMTIME_CHECK_START)/5.;
const double MEM_TIME = DT * LEN_TRAJ_K;//20.;
const int MEM_TIMESTEPS = MEM_TIME/DT; // number of memory timesteps
const double I_TIME = DT * LEN_TRAJ_I;//20.;
const int I_TIMESTEPS = I_TIME/DT; // number of memory timesteps
const double FINAL_TIME = MEM_TIME + DT; // final time of sigma(t)
const int FINAL_TIMESTEPS = FINAL_TIME/DT; // number of total timesteps
const int DOF_E_SQ_2 = 2 * DOF_E_SQ;
const int STAGES = 4; // stages, used within the VIDE algorithm
const int MAX_LOOPS_CONV = 20; // maximum number of loops for conv. algorithm
const double CONV_LIMIT = 1. * pow(10.,-2.); // value that the difference
// between sigma_z(t) and sigma_{z,max}(t) must be <= to be considered converged
const double CONV_ALG_BIG_STEP = 2.; // bigger time step within conv. alg.
const double CONV_ALG_SMALL_STEP = 0.25; // smaller time step within conv. alg.

#endif
