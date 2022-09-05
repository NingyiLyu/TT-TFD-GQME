# TT-TFD+GQME
This repository provides the code that carries out generalized quantum master equation (GQME) propagation, as well as the tensor-train thermo-field dynamics (TT-TFD) propagation to obtain the projection free input for the GQME routine. This code accompanies the paper "Tensor-Train Thermo-Field Memory Kernels for Generalized Quantum Master Equations" by Ningyi Lyu,* Ellen Mulvihill,* Micheline B. Soley, Eitan Geva, and Victor S. Batista, currently available on arxiv at [https://arxiv.org/abs/2208.14273v2](https://arxiv.org/abs/2208.14273v2).

# TT-TFD
The code in `SI_TT-TFD.py` generates the TT-TFD projection free inputs for spin-boson model 1 (detailed in section 5.1 of the main text). The current code initialize electronic state at donor state (line 70:tt_psi=tt_su), and the results correspond to four U elements (U_ab00) with a,b={0,1}. To obtain all 16 U elements, change line 65 to "tt_psi=tt_sd","tt_psi=tt_e1" and "tt_psi=tt_e2". Please refer to Section 4.4 of main text and Section S.IV of SI for details. With the current setting, each run should take a few hours to finish on a commercial computer. 

# Projection-Free Inputs (PFIs)

The projection-free input code in `PFIs_TT-TFD_GQME.py` will calculate the PFIs ${\cal F}(\tau)$ and $\dot{\cal F}(\tau)$ from the output of the TT-TFD code. It will calculate the linear combinations for the ${\cal U}(\tau)$ elements that start in the coherences [i.e., when $c \neq d$ in ${\cal U}_{abcd}(\tau)$], so there is no need for a separate code to do that. The code will run based on the model number and number of time steps specified at the beginning of the file and it will take the time step from the output of the TT-TFD code, so make sure those match what you just ran for the TT-TFD code. 

# Memory Kernels, Inhomogeneous Terms, and GQMEs

The code for the memory kernels, inhomogeneous terms, and GQMEs are all given in C++. They share a `constants.h` and `gen_constants.h` file, where you can specify the parameters of the model and GQME type. 

# Dependencies
 
 This package requires the following:
 
 - python
 - ttpy
 - matplotlib
 - numpy
 - scipy
 - c++11


 # Instructions
 
 ## TT-TFD Code
 
 The TT-TFD Code is given in `SI_TT-TFD.py`.
 
 To run the dynamics code, execute the following command:
 
 ```
 python SI_TT-TFD.py
 ```
 
### Output
The arrays "psu","psd","coh12" and "coh21" are the outputs for populations and coherences. For the current calculation (initialize electronic state at donor state, line 65:tt_psi=tt_su) should be directly used as U elements. For the ones that require linear combination, use the method described in SI, Section S.IV. 
 
## Projection-Free Inputs Code

The projection-free inputs (PFIs) code is given in `PFIs_TT-TFD_GQME.py`.

 To run the PFIs code, execute the following command:
 
 ```
 python PFIs_TT-TFD_GQME.py
 ```
## Memory Kernel Code

The memory kernel code is located in C++ > SubsetStates_GQME > KernelCode.

 To run the memory kernel code, after setting the gqmeType and number in the Makefile, execute the following commands:
 
 ```
 make
 ./<gqmeType>_main_Kernel<number>.exe
 ```
 
 ## Inhomogeneous Term Code

The inhomogeneous term code is located in C++ > SubsetStates_GQME > ITermCode.
 
  To run the inhomogeneous term code, after setting the gqmeType and number in the Makefile, execute the following commands:
 
 ```
 make
 ./<gqmeType>_main_ITerm<number>.exe
 ```
 
 ## GQME Code

The GQME code is located in C++ > SubsetStates_GQME > GQMECode.
 
To run the GQME code, after setting the gqmeType and number in the Makefile, execute the following commands:
 
 ```
 make
 ./<gqmeType>_main_GQME<number>.exe
 ```


