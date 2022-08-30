import numpy as np
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

#################
### Constants ###
#################
TIME_STEPS = 10 # number of time steps
MODEL_NUM = 1 # model number
DOF_E = 2 # number of electronic states
DOF_E_SQ = DOF_E * DOF_E

# setting parameters that are based on the model number
BETA = 5
GAMMA_DA = 1
if MODEL_NUM == 1:
    EPSILON = 1
    XI = 0.1
    OMEGA_C = 1
    OMEGA_MAX = 5
elif MODEL_NUM == 2:
    EPSILON = 1
    XI = 0.1
    OMEGA_C = 2
    OMEGA_MAX = 10
elif MODEL_NUM == 3:
    EPSILON = 1
    XI = 0.1
    OMEGA_C = 7.5
    OMEGA_MAX = 36
elif MODEL_NUM == 4:
    EPSILON = 1
    XI = 0.4
    OMEGA_C = 2
    OMEGA_MAX = 10
elif MODEL_NUM == 6:
    EPSILON = 0
    XI = 0.2
    OMEGA_C = 2.5
    OMEGA_MAX = 12

print("     time steps =", TIME_STEPS)
print("        model # =", MODEL_NUM)
print("          DOF_E =", DOF_E)
print("        epsilon =", EPSILON)
print("             xi =", XI)
print("        omega_c =", OMEGA_C)
print("      omega_max =", OMEGA_MAX)

##############################################################
### Sections Involving the Calculation of ${\cal U}(\tau)$ ###
##############################################################

INIT_NAMES = ["initsu", "inite1", "inite2", "initsd"]
INIT_NUMS = ["00", "01sx", "10sy", "11"]
STATE_NAMES = ["psu", "coh21", "coh12", "psd"]
STATE_NUMS = ["00", "01", "10", "11"]
FILE_PREFIX = "%s_steps_"%TIME_STEPS

# matrix for U
U = np.zeros((TIME_STEPS, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)

# Converting TT-TFD Files to ${\cal U}(\tau)$ Matrix
initNum = 0 # number that is advanced in each loop for the init state
for initName in INIT_NAMES:
    # pull in time steps
    infileStrTime = "Output/TT-TFD_Output/t.npy"
    time = np.load(infileStrTime)
    
    stateNum = 0 # number that is advanced in each loop for the dynamic state
    for stateName in STATE_NAMES:
        # initializing the real and imag parts of U
        real = np.zeros((TIME_STEPS))
        imag = np.zeros((TIME_STEPS))
        
        # string to load in dynamics from TT-TFD
        infileStr = "Output/TT-TFD_Output/" + stateName + "_" + initName + ".npy"
        
        # if dynamic state is a population, files only have real part
        if stateNum == 0 or stateNum == 3:
            real = np.load(infileStr)
            
        # if dynamic state is a coherence, files have real and imag part
        else:
            numbers = np.load(infileStr)
            
            for i in range(0,TIME_STEPS):
                real[i] = np.real(numbers[i])
                imag[i] = np.imag(numbers[i])
        
        for i in range(TIME_STEPS):
            U[i][stateNum][initNum] = real[i] + 1.j * imag[i]
        
        stateNum += 1
    initNum += 1

# Calculating the Coherence Elements of ${\cal U}(\tau)$
for j in range(DOF_E_SQ):
    a = str(int(j/DOF_E)) # a index of U_{abcd}
    b = str(int(j%DOF_E)) # b index of U_{abcd}
    for cc in range(DOF_E):
        c = str(cc) # c index of U_{abcd}
        for dd in range(cc + 1, DOF_E):
            d = str(dd) # d index of U_{abcd}
            # creates the matching vector index from c and d
            k = DOF_E * cc + dd
            # creates the matching vector index from d and c
            oppositeIndex = DOF_E * dd + cc
            # creates the matching vector index from c and c
            DDIndex = DOF_E * cc + cc
            # creates the matching vector index from d and d
            AAIndex = DOF_E * dd + dd
                
            for i in range(TIME_STEPS):
                # holds the abcd_sx value needed for calculating abdc that is
                # overwritten when calculating abcd first
                hold = U[i][j][k]
                
                # calculates abcd = abcd_sx + i * abdc_sy
                #                         - 1/2 * (1 + i) * (abcc + abdd)
                U[i][j][k] = U[i][j][k] + 1.j * U[i][j][oppositeIndex]
                U[i][j][k] -= 0.5 * (1. + 1.j) * (U[i][j][DDIndex] + U[i][j][AAIndex])
                
                # calculates abdc = abcd_sx - i * abdc_sy
                #                         - 1/2 * (1 - i) * (abcc + abdd)
                U[i][j][oppositeIndex] = hold - 1.j * U[i][j][oppositeIndex]
                U[i][j][oppositeIndex] -= 0.5 * (1. - 1.j) * (U[i][j][DDIndex] + U[i][j][AAIndex])


# Printing ${\cal U}(\tau)$ Matrix to Files

for j in range(DOF_E_SQ):
    a = str(int(j/DOF_E))
    b = str(int(j%DOF_E))
    for k in range(DOF_E_SQ):
        c = str(int(k/DOF_E))
        d = str(int(k%DOF_E))

        # outfile for the U
        outfileStr = "Output/U_Output/U_" + a + b + c + d + "_TT-TFD_"
        outfileStr += FILE_PREFIX + "model%s"%MODEL_NUM  + ".dat"
        #print("\t   ", outfileStr)
        f = open(outfileStr, "w")

        for i in range(TIME_STEPS):
            f.write("%s\t%s\t%s\n"%(time[i], U[i][j][k].real, U[i][j][k].imag))
        f.close()

#############################################################
### Calculating ${\cal F}(\tau)$ and $\dot{\cal F}(\tau)$ ###
#############################################################

# variables
time = np.zeros((TIME_STEPS))
U = np.zeros((TIME_STEPS, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)
F = np.zeros((TIME_STEPS, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)
Fdot = np.zeros((TIME_STEPS, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)

# Calling Time and ${\cal U}(\tau)$ Values From Files
for j in range(DOF_E_SQ):
    a = str(int(j/DOF_E))
    b = str(int(j%DOF_E))
    for k in range(DOF_E_SQ):
        c = str(int(k/DOF_E))
        d = str(int(k%DOF_E))

        t, Ureal, Uimag = np.hsplit(
                np.loadtxt("Output/U_Output/U_%s%s%s%s_TT-TFD_"%(a,b,c,d) 
                + "%s_steps_model%s.dat"%(TIME_STEPS, MODEL_NUM)), 3)

        for i in range(TIME_STEPS):
            time[i] = t[i]
            U[i][j][k] = Ureal[i] + 1.j * Uimag[i]

# constants that depend on the time steps
DT = float(time[1]) # time step
FINAL_TIME = float(time[len(time) - 1]) # final time

# setting parameter string
PARAM_STR = "_Ohmic_TT-TFD_b5G1_e%s_t%.8f_"%(EPSILON, DT)
PARAM_STR += "xi%swc%s_wmax%s_dofn60_tf%.4f"%(XI, OMEGA_C, OMEGA_MAX, FINAL_TIME)

print("             DT =", DT)
print("     final time =", FINAL_TIME)
print("   param string =", PARAM_STR)

# Function to Print ${\cal F}(\tau)$ and $\dot{\cal F}(\tau)$
def printFFdot(abcdStr, timeFFdot, Freal, Fimag, Fdotreal, Fdotimag):
    # outfileStr for F
    outfileFStr = "Output/ProjFree_Output/F_"
    outfileFStr += abcdStr + PARAM_STR + ".dat"
    
    f = open(outfileFStr, "w")
    
    for i in range(len(Freal)):
        f.write("%s\t%s\t%s\n"%(timeFFdot[i], Freal[i], Fimag[i]))
    f.close()
    
    # outfileStr for Fdot
    outfileFdotStr = "Output/ProjFree_Output/Fdot_"
    outfileFdotStr += abcdStr + PARAM_STR + ".dat"
                
    f = open(outfileFdotStr, "w")
    for i in range(len(Fdotreal)):
        f.write("%s\t%s\t%s\n"%(timeFFdot[i], Fdotreal[i], Fdotimag[i]))
    f.close()


# 2nd-Order Central Difference
for j in range(DOF_E_SQ):
    a = str(int(j/DOF_E)) # a index of PFI_{abcd}
    b = str(int(j%DOF_E)) # b index of PFI_{abcd}
    for k in range(DOF_E_SQ):
        c = str(int(k/DOF_E)) # c index of PFI_{abcd}
        d = str(int(k%DOF_E)) # d index of PFI_{abcd}
            
        # derivative of U for F, F = i\dot{U} so
        # Freal = -1 * \dot{Uimag} and Fimag = \dot{Ureal}
        Ureal = np.zeros((TIME_STEPS))
        Uimag = np.zeros((TIME_STEPS))
        for i in range(TIME_STEPS):
            Ureal[i] = U[i][j][k].real
            Uimag[i] = U[i][j][k].imag
                
        Freal = -1. * np.gradient(Uimag.flatten(), DT, edge_order = 2)
        Fimag = np.gradient(Ureal.flatten(), DT, edge_order = 2)
                
        # first derivative of F for Fdot, Fdot = \dot{F} so
        # Fdotreal = \dot{Freal} and Fdotimag = \dot{Fimag}
        Fdotreal = np.gradient(Freal, DT)
        Fdotimag = np.gradient(Fimag, DT)
                
        abcdStr = a + b + c + d
        printFFdot(abcdStr, time, Freal, Fimag, Fdotreal, Fdotimag)
