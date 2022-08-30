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
#include "FunctionTemplates_Kernel.h"

// Purpose: For the first function, the purpose is to calculate an integral of
//          the form:
//            kernel = linearTerm + i * \int_0^tau dtau' F(tau - tau')
//                                                           * kernel(tau')
//          involving Complex_3D_Matrix via trapezoid method
//
//          For the second function, the purpose is to run the Volterra
//          convergence algorithm and print the converged kernel
//
// Created by Ellen Mulvihill on 8/29/22
// Based on code created by Ellen Mulvihill on 9/26/18
//
// Most recent modification on 8/29/22 (Ellen)

// calculates an integral of the form
// kernel(t_f) = linearTerm(t_f)
//                  + i*\int_(t_0)^(t_f) dt' F(t_f - t') * kernel(t')
void CalculateIntegral(int& numStates, Complex_3D_Matrix& linearTerm,
                       Complex_3D_Matrix& F, Complex_3D_Matrix& prevKernel,
                       Complex_3D_Matrix& kernel){
    int i,j,k,l,n,c; // integers for loops
    
    // loops to calculate kernel, with integral calculated via extended
    // trapezoidal method then multiplied by i and added to the linearTerm
    for (n = 1; n < LEN_TRAJ_K; n++){ // t_0 already calculated, so starts at 1
        for (i = 0; i < numStates; i++){
            for (j = 0; j < numStates; j++){
                kernel[n][i][j] = 0.;
                for (k = 0; k < numStates; k++){
                    kernel[n][i][j] += 0.5 * INTEGRAL_STEP_K * F[n][i][k] * kernel[0][k][j];
                    kernel[n][i][j] += 0.5 * INTEGRAL_STEP_K * F[0][i][k] * prevKernel[n][k][j];
                    for (c = 1; c < n; c++){
                        // since a new (supposed-to-be-better) guess for the
                        // kernel has been calculated for previous time
                        // steps, can use it rather than prevKernel
                        kernel[n][i][j] += INTEGRAL_STEP_K * F[n-c][i][k] * kernel[c][k][j];
                    }
                }
                kernel[n][i][j] = I * kernel[n][i][j] + linearTerm[n][i][j];
            }
        }
    }
    
    return;
}

// Calculates a Volterra equation of the form
// kernel(t_f) = linearTerm(t_f)
//                + i*\int_(t_0)^(t_f) dt' integralInput(t_f - t') * kernel(t')
// iteratively with a set convergence parameter and maximum number of
// iterations, given in ../../constants.h
void RunVolterra(string& idStr_K, vector<string>& states, int& numStates,
                 Complex_3D_Matrix& linearTerm, Complex_3D_Matrix& F){
    int i,j,k,l,n; // integers for loops
    int numIter; // integer for the number of the interation
    string outputStr("");
    ofstream outfile; // variable for output files
    
    // start timing for iterations
    clock_t clock_start;
    clock_start = clock();
    
    // matrix holding the previous guess of the kernel
    Complex_3D_Matrix prevKernel(LEN_TRAJ_K, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
    Complex_3D_Matrix kernel(LEN_TRAJ_K, Complex_Matrix(numStates, vector<Complex>(numStates, 0.0)));
    
    // sets first guess of the kernel to the linear term (aka everything not in
    // the integral)
    for (i = 0; i < numStates; i++){
        for (j = 0; j < numStates; j++){
            for (n = 0; n < LEN_TRAJ_K; n++){
                prevKernel[n][i][j] = linearTerm[n][i][j];
            }
        }
    }
    
    // sets the kernel(t_0) = linearTerm(t_0) (since integral is zero at t_0)
    for (i = 0; i < numStates; i++){
        for (j = 0; j < numStates; j++){
            kernel[0][i][j] = linearTerm[0][i][j];
        }
    }
    
    // opens the file to print the info about the volterra iterations.
    // If the iterations succeed, will print the number of iterations to
    // converge the kernel and the number of iterations to converge all the
    // error bars and the time it took total.
    // If the iterations fail, will print a fail message; the prevKernel and
    // kernel values for each time step, abcd indices, and error bars; and
    // the time it took
    outputStr = "VolterraInfo_K_" + GQME_TYPE + "_";
    //outputStr = "VolterraInfo_K_blockBLOckNUM_" + GQME_TYPE + "_";
    if (GQME_TYPE == "SubsetStates" || GQME_TYPE == "SingleState"){
        for (i = 0; i < numStates; i++){
            outputStr += states[i] + "_";
        }
    }
    outputStr += BATH_TYPE + "_" + DYN_STRING + "_";
    
    outfile.open((K_FOLDER + outputStr + idStr_K + ".dat").c_str());
    if (!outfile.is_open()) { // prints error if file not open
        cout << "\t ERROR: output file of K " << GQME_TYPE;
        cout << " Volterra information cannot open" << endl;
        cout << K_FOLDER + outputStr + idStr_K + ".dat" << endl;
    }
    cout << "\t Convergence Parameter: " << CONVERGENCE_PARAM << endl;
    outfile << "\t Convergence Parameter: " << CONVERGENCE_PARAM << endl;
    
    // starts iterating over numIter up to MAX_ITERS specified in ../constants.h
    for (numIter = 0; numIter < MAX_ITERS; numIter++){
        // prints the iteration number
        cout << "\t Iteration: " << numIter << endl;
        outfile << "\t Iteration: " << numIter << endl;

        // calculates the Volterra equation for this iteration, with the
        // function given above in this file
        CalculateIntegral(numStates, linearTerm, F, prevKernel, kernel);
        
        // vector of convergence for the kernel and each error bar
        int numConv(0);
        for (i = 0; i < numStates; i++){
            for (j = 0; j < numStates; j++){
                for (n = 0; n < LEN_TRAJ_K; n++){
                    // For each time step and matrix element, if the kernel and
                    // prevKernel are  within the convergence parameter set in
                    // ../../constants.h, then 1 is added to numConv
                    if (abs(kernel[n][i][j] - prevKernel[n][i][j]) <= CONVERGENCE_PARAM){
                        numConv += 1;
                    }
                    else {
                        // if the maximum number of iterations has been reached,
                        // prints the kernel and prevKernel values for the time
                        // steps and matrix elements that did not converge at
                        // the maximum iteration to VolterraInfo
                        if (numIter == MAX_ITERS - 1){
                            outfile << "\t\t K time step and matrix element ";
                            outfile << "that didn't converge: ";
                            outfile << n << ", " << i << j << endl;
                            outfile << "\t\t prevKernel: ";
                            outfile << prevKernel[n][i][j] << endl;
                            outfile << "\t\t          K: ";
                            outfile << kernel[n][i][j] << endl;
                        }
                    }
                }
            }
        }
        // For each time step, if it is converged, 1 is added. Therefore, the
        // kernel is converged when numConv = LEN_TRAJ
        if (numConv == LEN_TRAJ_K * numStates * numStates){
            // prints the number of iterations to the shell and VolterraInfo
            cout << "\t >>> Number of iterations: " << numIter << endl;
            outfile << "Number of iterations: " << numIter << endl;
                
            // outputs the converged kernel. Function is located in
            // PrintKernel.cpp
            PrintKernel("", idStr_K, states, numStates, kernel);
            
            //cout << "linearTerm(0) = " << linearTerm[0][0][0] << endl;
            //cout << "kernel(0) = " << kernel[0][0][0] << endl;
                
            double ptc, pts, ptm, pth;
            ptc = clock() - clock_start;
            pts = ptc/CLOCKS_PER_SEC;
            ptm = pts/60.;
            pth = ptm/60.;
            if (ptm < 2.){ // if total time < 2 minutes, prints seconds
                cout << "\t Volterra time: " << pts << " seconds" << endl;
                outfile << "Volterra time: " << pts << " seconds" << endl;
            }
            else if (pth < 2.){ // if total time < 2 hours, prints minutes
                cout << "\t Volterra time: " << ptm << " minutes" << endl;
                outfile << "Volterra time: " << ptm << " minutes" << endl;
            }
            else { // if total time >= 2 hours, prints hours
                cout << "\t Volterra time: " << pth << " hours" << endl;
                outfile << "Volterra time: " << pth << " hours" << endl;
            }
                
            // closes the VolterraInfo file and frees up outfile
            outfile.close();
            outfile.clear();
                
            // exits the RunVolterra function since everything has converged
            break;
        }
        
        // Before starting a new iteration, sets prevKernel to kernel and
        // kernel to 0
        for (i = 0; i < numStates; i++){
            for (j = 0; j < numStates; j++){
                for (n = 1; n < LEN_TRAJ_K; n++){ // kernel is defined at t = 0
                    prevKernel[n][i][j] = kernel[n][i][j];
                    kernel[n][i][j] = 0.;
                }
            }
        }
        
        // prints the time each iteration took
        if (((clock()-clock_start)/(double) CLOCKS_PER_SEC)/60. < 120){
            printf("\t Iteration runtime: %5.8f minutes.\n",((clock()-clock_start)/(double) CLOCKS_PER_SEC)/60.);
            outfile << "\t Iteration runtime: " << ((clock()-clock_start)/(double) CLOCKS_PER_SEC)/60.;
            outfile << " minutes." << endl;
        }
        else {
            printf("\t Iteration runtime: %5.8f hours.\n",((clock()-clock_start)/(double) CLOCKS_PER_SEC)/60./60.);
            outfile << "\t Iteration runtime: "<< ((clock()-clock_start)/(double) CLOCKS_PER_SEC)/60./60.;
            outfile << " hours." << endl;
        }
        
        // if the maximum number of iterations has been reached, an error
        // message is printed to the shell and VolterraInfo file
        if (numIter == MAX_ITERS - 1){
            cout << "\t ERROR: Did not converge for " << MAX_ITERS;
            cout << " iterations" << endl;
            outfile << "ERROR: Did not converge for " << MAX_ITERS;
            outfile << " iterations" << endl;
            
            // prints time of the RunVolterra function to the shell and
            // VolterraInfo file
            double ptc, pts, ptm, pth;
            ptc = clock() - clock_start;
            pts = ptc/CLOCKS_PER_SEC;
            ptm = pts/60.;
            pth = ptm/60.;
            if (ptm < 2.){ // if total time < 2 minutes, prints seconds
                cout << "\t Volterra time: " << pts << " seconds" << endl;
                outfile << "Volterra time: " << pts << " seconds" << endl;
            }
            else if (pth < 2.){ // if total time < 2 hours, prints minutes
                cout << "\t Volterra time: " << ptm << " minutes" << endl;
                outfile << "Volterra time: " << ptm << " minutes" << endl;
            }
            else { // if total time >= 2 hours, prints hours
                cout << "\t Volterra time: " << pth << " hours" << endl;
                outfile << "Volterra time: " << pth << " hours" << endl;
            }
            
            // closes the VolterraInfo file and clears up outfile
            outfile.close();
            outfile.clear();
        }
    }
    
    return;
}
