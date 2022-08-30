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
TIME_STEPS = 10200 # number of time steps
MODEL_NUM = 5 # model number
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
elif MODEL_NUM == 5:
    BETA = 3
    GAMMA_DA = "0.333333"
    EPSILON = 0
    XI = 0.1
    OMEGA_C = 1
    OMEGA_MAX = 5
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

############################
### Functions for Graphs ###
############################

## 4x4 Graph for ${\cal U}(\tau)$, ${\cal F}(\tau)$, or $\dot{\cal F}(\tau)$
def graph4x4(real_imag, modelNum, methods, numCols, legendPos_x, legendPos_y, 
             quantityStr, graphTitleStr, graphStr, timeDict, quantityDict):
    linestyles = ["-", "--", "-.", ":", ":", ":", "-.", ":"]
    colors = ['b', 'r', 'm', 'g', 'c', 'y']
    linewidths = [2,2,2,2,3,3]

    fig = plt.figure(figsize = (18,10.5))

    # creates the 4x4 graphs
    ax = [] 
    for i in range(0,4):
        for j in range(0,4):
            ax.append(plt.subplot2grid((4, 4), (i, j)))

    # sets the spacing between plots
    plt.subplots_adjust(wspace = 0.85, hspace = 0)
  
    for k in range(0, methods):
        # pulling the time and quantity from the dictionaries 
        time = timeDict[graphStr[k]]
        quantity = quantityDict[graphStr[k]]

        # look at real or imag part depending on real_imag input
        if real_imag == "Real":
            quantity = quantity.real
        elif real_imag == "Imag":
            quantity = quantity.imag
        else:
            print("ERROR: real_imag value not Real or Imag")
            break

        # making sure the time values cut off at the limit of the quantity
        # so that their lengths match for the plot
        rangeLimit = len(quantity[:,0,0])  

        # loops to plot graphs
        for j in range(0, 4):
            # ab indices of quantity_{abcd}
            l = str(int(j/DOF_E)) + str(int(j%DOF_E)) 

            # plotting the quantity
            ax[0 + j*4].plot(time[0:rangeLimit], quantity[:,j,0], 
                             color = colors[k], linewidth=linewidths[k], 
                             linestyle=linestyles[k], label = r'%s'%graphStr[k])
            ax[1 + j*4].plot(time[0:rangeLimit], quantity[:,j,1], 
                             color = colors[k], linewidth=linewidths[k], 
                             linestyle=linestyles[k])
            ax[2 + j*4].plot(time[0:rangeLimit], quantity[:,j,2], 
                             color = colors[k], linewidth=linewidths[k], 
                             linestyle=linestyles[k])
            ax[3 + j*4].plot(time[0:rangeLimit], quantity[:,j,3], 
                             color = colors[k], linewidth=linewidths[k], 
                             linestyle=linestyles[k])
            
            # top 3 rows of graphs
            if j < 3: 
                # since the graphs share x-axes, we need to turn off the ticks 
                # for the upper three graphs in each column
                ax[0 + j*4].set(xticks=[])
                ax[1 + j*4].set(xticks=[])
                ax[2 + j*4].set(xticks=[])
                ax[3 + j*4].set(xticks=[])

                # makes the y tick values larger
                ax[0 + j*4].tick_params(axis='y', labelsize=16)
                ax[1 + j*4].tick_params(axis='y', labelsize=16)
                ax[2 + j*4].tick_params(axis='y', labelsize=16)
                ax[3 + j*4].tick_params(axis='y', labelsize=16)

                # controls the number of y ticks
                ax[0 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
                ax[1 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
                ax[2 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
                ax[3 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
                
            # bottom row of graphs
            else: 
                # makes both tick values larger 
                ax[0 + j*4].tick_params(axis='both', labelsize=16)
                ax[1 + j*4].tick_params(axis='both', labelsize=16)
                ax[2 + j*4].tick_params(axis='both', labelsize=16)
                ax[3 + j*4].tick_params(axis='both', labelsize=16)

                # controls the number of ticks
                ax[0 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[1 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[2 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[3 + j*4].yaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[0 + j*4].xaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[1 + j*4].xaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[2 + j*4].xaxis.set_major_locator(MaxNLocator(nbins=5))
                ax[3 + j*4].xaxis.set_major_locator(MaxNLocator(nbins=5))

            # y labels for U and F
            if len(quantityStr) == 1:
                ax[0 + j*4].set_ylabel(r'${\cal %s}_{%s00}$'%(quantityStr, l),fontsize = 28)
                ax[1 + j*4].set_ylabel(r'${\cal %s}_{%s01}$'%(quantityStr, l),fontsize = 28)
                ax[2 + j*4].set_ylabel(r'${\cal %s}_{%s10}$'%(quantityStr, l),fontsize = 28)
                ax[3 + j*4].set_ylabel(r'${\cal %s}_{%s11}$'%(quantityStr, l),fontsize = 28)

            # y labels for Fdot
            elif quantityStr == "Fdot": 
                ax[0 + j*4].set_ylabel(r'$\dot{\cal F}_{%s00}$'%(l),fontsize = 28)
                ax[1 + j*4].set_ylabel(r'$\dot{\cal F}_{%s01}$'%(l),fontsize = 28)
                ax[2 + j*4].set_ylabel(r'$\dot{\cal F}_{%s10}$'%(l),fontsize = 28)
                ax[3 + j*4].set_ylabel(r'$\dot{\cal F}_{%s11}$'%(l),fontsize = 28)

            elif quantityStr[2] == "i":
                ax[0 + j*4].set_ylabel(r'${\cal %s}^{diff}_{%s00}$'%(quantityStr[0], l),fontsize = 28)
                ax[1 + j*4].set_ylabel(r'${\cal %s}^{diff}_{%s01}$'%(quantityStr[0], l),fontsize = 28)
                ax[2 + j*4].set_ylabel(r'${\cal %s}^{diff}_{%s10}$'%(quantityStr[0], l),fontsize = 28)
                ax[3 + j*4].set_ylabel(r'${\cal %s}^{diff}_{%s11}$'%(quantityStr[0], l),fontsize = 28)

            elif quantityStr[4] == "d":
                ax[0 + j*4].set_ylabel(r'$\dot{\cal F}^{diff}_{%s00}$'%(l),fontsize = 28)
                ax[1 + j*4].set_ylabel(r'$\dot{\cal F}^{diff}_{%s01}$'%(l),fontsize = 28)
                ax[2 + j*4].set_ylabel(r'$\dot{\cal F}^{diff}_{%s10}$'%(l),fontsize = 28)
                ax[3 + j*4].set_ylabel(r'$\dot{\cal F}^{diff}_{%s11}$'%(l),fontsize = 28)
                

    # sets x labels on bottom row of graphs
    ax[12].set_xlabel(r'$\Gamma\tau$',fontsize = 32)
    ax[13].set_xlabel(r'$\Gamma\tau$',fontsize = 32)
    ax[14].set_xlabel(r'$\Gamma\tau$',fontsize = 32)
    ax[15].set_xlabel(r'$\Gamma\tau$',fontsize = 32)
    
    # puts a buffer on the left side, as y labels have been cut off before 
    plt.gcf().subplots_adjust(left=0.15)

    # gets the labels from the first graph
    handles,labels = ax[0].get_legend_handles_labels()
    legend = ax[0].legend(handles,labels,loc = 'upper right', 
                          bbox_to_anchor=(legendPos_x * 4.6, legendPos_y * 1.5), 
                          fontsize = 20, borderpad=0.2, borderaxespad=0.2, 
                          ncol = numCols, 
                          title = "%s, %s"%(real_imag, graphTitleStr))
    # adjusts title settings 
    plt.setp(legend.get_title(),fontsize = 20, fontweight='bold')

    # generates a string to differentiate the files
    outputStr = "_%s_"%(real_imag) + graphTitleStr
    for k in range(methods):
        outputStr += "_" + graphStr[k] 
    
    # saves the figure
    plt.savefig("Figures/" + quantityStr + "_Graphs/" + quantityStr 
                + "_model%s"%(modelNum) + outputStr + ".pdf", dpi=plt.gcf().dpi, 
                bbox_inches='tight')

##############################################################
### Sections Involving the Calculation of ${\cal U}(\tau)$ ###
##############################################################

INIT_NAMES = ["initsu", "initcoh1", "initcoh2", "initsd"]
INIT_NUMS = ["00", "01sx", "10sy", "11"]
STATE_NAMES = ["psu", "cohdu", "cohud", "psd"]
STATE_NUMS = ["00", "01", "10", "11"]
FILE_PREFIX = "%s_steps_"%TIME_STEPS

# matrix for U
U = np.zeros((TIME_STEPS, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)

# Converting TT-TFD Files to ${\cal U}(\tau)$ Matrix
initNum = 0 # number that is advanced in each loop for the init state
for initName in INIT_NAMES:
    # pull in time steps
    infileStrTime = "Output/TT-TFD_Output/" + FILE_PREFIX + "model%s"%MODEL_NUM + "_" 
    infileStrTime += initName + "/" + FILE_PREFIX + "model%s"%MODEL_NUM  + "_" 
    infileStrTime += initName + "_time.csv"
    time = np.loadtxt(infileStrTime, dtype=np.dtype('f8'))
    
    stateNum = 0 # number that is advanced in each loop for the dynamic state
    for stateName in STATE_NAMES:
        # initializing the real and imag parts of U
        real = np.zeros((TIME_STEPS))
        imag = np.zeros((TIME_STEPS))
        
        # string to load in dynamics from TT-TFD
        infileStr = "TT-TFD_Output/" + FILE_PREFIX + "model%s"%MODEL_NUM + "_" 
        infileStr += initName + "/" + FILE_PREFIX + "model%s"%MODEL_NUM + "_" 
        infileStr += initName + "_" + stateName + ".csv"
        
        # if dynamic state is a population, files only have real part
        if stateNum == 0 or stateNum == 3:
            real = np.loadtxt(infileStr, dtype=np.dtype('f8'))
            
        # if dynamic state is a coherence, files have real and imag part
        else:
            numbers = np.loadtxt(infileStr, dtype=np.complex_)
            
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
        outfileStr = "U_Output/U_" + a + b + c + d + "_TT-TFD_"
        outfileStr += FILE_PREFIX + "model%s"%MODEL_NUM  + ".dat"
        #print("\t   ", outfileStr)
        f = open(outfileStr, "w")

        for i in range(TIME_STEPS):
            f.write("%s\t%s\t%s\n"%(time[i], U[i][j][k].real, U[i][j][k].imag))
        f.close()

# Graphing ${\cal U}(\tau)$ from Files

MODEL_NUM = 5
timeStepOptions = [10200]#[2000, 3400, 5100, 10000]
graphStr = ["10200"]#["2000", "3400", "5100", "10000"]
methods = 1
numCols = 1

UTimeDict_file = {}
UDict_file = {}

for timeSteps in timeStepOptions:
    U_file = np.zeros((timeSteps, DOF_E_SQ, DOF_E_SQ), dtype=np.complex_)
    time_file = np.zeros((timeSteps))
    for j in range(DOF_E_SQ):
        a = int(j/DOF_E)
        b = int(j%DOF_E)
        for k in range(DOF_E_SQ):
            c = int(k/DOF_E)
            d = int(k%DOF_E)

            t, Ureal, Uimag = np.hsplit(
                np.loadtxt("U_Output/U_%s%s%s%s_TT-TFD_"%(a,b,c,d) 
                + "%s_steps_model%s.dat"%(timeSteps, MODEL_NUM)), 3)

            for i in range(timeSteps):
                time_file[i] = t[i]
                U_file[i][j][k] = Ureal[i] + 1.j * Uimag[i]
    
    UTimeDict_file.update({"%s"%timeSteps : time_file})
    UDict_file.update({"%s"%timeSteps : U_file})


graphTitleStr = "Model %s"%MODEL_NUM#"Chatterjee & Makri"
legendPos_x = 1#0.75
legendPos_y = 1.
graph4x4("Real", MODEL_NUM, methods, numCols, legendPos_x, legendPos_y, "U", 
         graphTitleStr, graphStr, UTimeDict_file, UDict_file)
graph4x4("Imag", MODEL_NUM, methods, numCols, legendPos_x, legendPos_y, "U", 
         graphTitleStr, graphStr, UTimeDict_file, UDict_file)

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

        if CandM_BOOL == True:
            t, Ureal, Uimag = np.hsplit(
                np.loadtxt("U_Output/U_%s%s%s%s_TT-TFD_"%(a,b,c,d) 
                + "%s_steps_makri_model%s.dat"%(TIME_STEPS, MODEL_NUM)), 3)
        else:
            t, Ureal, Uimag = np.hsplit(
                np.loadtxt("U_Output/U_%s%s%s%s_TT-TFD_"%(a,b,c,d) 
                + "%s_steps_model%s.dat"%(TIME_STEPS, MODEL_NUM)), 3)

        for i in range(TIME_STEPS):
            time[i] = t[i]
            U[i][j][k] = Ureal[i] + 1.j * Uimag[i]

# constants that depend on the time steps
DT = float(time[1]) # time step
FINAL_TIME = float(time[len(time) - 1]) # final time

# setting parameter string
if TIME_STEPS == 3400 or TIME_STEPS == 10200:
    PARAM_STR = "_Ohmic_TT-TFD_b%sG%s_e%s_t%.8f_"%(BETA, GAMMA_DA, EPSILON, DT)
    PARAM_STR += "xi%swc%s_wmax%s_dofn60_tf%.3f"%(XI, OMEGA_C, OMEGA_MAX, FINAL_TIME)
else:
    PARAM_STR = "_Ohmic_TT-TFD_b5G1_e%s_t%.8f_"%(EPSILON, DT)
    PARAM_STR += "xi%swc%s_wmax%s_dofn60_tf%.4f"%(XI, OMEGA_C, OMEGA_MAX, FINAL_TIME)

print("             DT =", DT)
print("     final time =", FINAL_TIME)
print("   param string =", PARAM_STR)

# Function to Print ${\cal F}(\tau)$ and $\dot{\cal F}(\tau)$
def printFFdot(abcdStr, timeFFdot, Freal, Fimag, Fdotreal, Fdotimag):
    # outfileStr for F
    outfileFStr = "ProjFree_Output/F_"
    outfileFStr += abcdStr + PARAM_STR + ".dat"
    
    f = open(outfileFStr, "w")
    
    for i in range(len(Freal)):
        f.write("%s\t%s\t%s\n"%(timeFFdot[i], Freal[i], Fimag[i]))
    f.close()
    
    # outfileStr for Fdot
    outfileFdotStr = "ProjFree_Output/Fdot_"
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
