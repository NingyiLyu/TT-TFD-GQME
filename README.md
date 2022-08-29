# TT-TFD_GQME
This repository provides the code that carries out generalized quantum master equation (GQME) propagation, as well as the tensor-train thermo-field dynamics (TT-TFD) propagation to obtain the projection free input for the GQME routine. 

# TT-TFD
The code in SI_TT-TFD.py generates the TT-TFD projection free inputs for spin-boson model 1 (detailed in section 5.1 of the main text). These would correspond to four U elements (U_ab00) with a,b={0,1}. To obtain all 16 U elements, change line 70 to "tt_psi=tt_sd","tt_psi=tt_e1" and "tt_psi=tt_e2". Please refer to Section 4.4 of main text and Section S.IV of SI for details.

# Dependencies
 
 This package requires the following:
 
 - python
 - ttpy
 - matplotlib
 - numpy
 - scipy

 
 # Instructions
 
 To run the dynamics code, execute the following commands:
 
 ```
 python SI_TT-TFD.py
 ```

# Output
The arrays "psu","psd","coh12" and "coh21" are the outputs for populations and coherences. For the current calculation (initialize electronic state at donor state, line 70:tt_psi=tt_su) should be directly used as U elements. For the ones that require linear combination, use the method described in SI, Section S.IV. 
