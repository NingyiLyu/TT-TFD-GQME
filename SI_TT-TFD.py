#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 04:11:14 2021

@author: ningyi
"""
import numpy as np
import tt
import tt.ksl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from time import process_time
t_start=process_time()
qmodes=60            #number of quantum bath modes
nsc = 5100             # number of propagation steps
tau = 10              # propagation time step
eps = 1e-14            # tt approx error
rma = 2000000                # max tt rank
dim = qmodes         # number of coords
nstates=2              # number of surfaces
occ=10                 #maximum occupation number; low for harmonic systems
EYE=1j
kondo=0.1            #kondo parameter
au=3.00166*10**(-4)
cmn1toau=4.5563353e-6 # Conversion of wavenumbers to atomic units
au2ps=0.00002418884254# Conversion of attoseconds to atomic units
wc=1.*au     #max freq for Ohmic bath discretization
wmax=5.*au
om=wc/qmodes*(1-np.exp(-wmax/wc))
J2au=2.2937104486906*10**17 #Conversion of Joules to atomic units
kB=1.308649*10**(-23) #Boltzmann's constant, in J/K
T=300                #Temperature, in K
kT=0.2*au         #kBT, in atomic units
beta=1./kT            #beta
lam=20.*cmn1toau
alpha=lam/(2*wc*np.pi)
#initialize arrays for parameters
freq=np.zeros((qmodes)) #frequency
ck=np.zeros((qmodes))   #linear electron-phonon coupling constant
gk=np.zeros((qmodes))   #ck in occupation number representation
thetak=np.zeros((qmodes)) #temperature-dependent mixing parameter in TFD
sinhthetak=np.zeros((qmodes)) #sinh(theta)
coshthetak=np.zeros((qmodes)) #cosh(theta)
for i in range(qmodes):                   
    freq[i]=-wc*np.log(1-(i+1)*om/(wc)) # Ohmic frequency
    ck[i]=np.sqrt(kondo*om)*freq[i] #Ohmic coupling constant
    gk[i]=-ck[i]/np.sqrt(2*freq[i]) #Transfer ck to occ. num. representation
    thetak[i]=np.arctanh(np.exp(-beta*freq[i]/2)) #theta, defined for harmonic models
    sinhthetak[i]=np.sinh(thetak[i]) #sinh(theta)
    coshthetak[i]=np.cosh(thetak[i]) #cosh(theta)
eelec=1.0*au #electronic state energy, in a.u.
coupling=1.0*au #electronic interstate coupling, in a.u.
#Build initial ground state
su=np.array([1,0]) 
sd=np.array([0,1])
e1=np.sqrt(0.5)*(su+sd)
e2=np.sqrt(0.5)*(su+EYE*sd)
tt_su=tt.tensor(su)
tt_sd=tt.tensor(sd)
tt_e1=tt.tensor(e1)
tt_e2=tt.tensor(e2)
tt_Ie=tt.eye(2,1)
gs=np.zeros((occ))
gs[0]=1.
tt_gs=tt.tensor(gs)
tt_psi0=tt_su
for k in range(2*qmodes):#double space formation
    tt_psi0=tt.kron(tt_psi0,tt_gs)
#constructing Pauli operators
px=np.array([[0,1],[1,0]])
pz=np.array([[1,0],[0,-1]])
#Build electronic site energy matrix
He=eelec*pz-coupling*px
#TT-ize that energy matrix
tt_He=tt.matrix(He)
tt_He=tt.kron(tt_He,tt.eye(occ,qmodes*2))
#Build number operator, corresponds to harmonic oscillator Hamiltonian, see notesc
numoc=np.diag(np.arange(0,occ,1))
tt_numoc=tt.eye(occ,qmodes)*0.
#Build displacement operator, corresponds to x operator in real space
eneroc=np.zeros((occ,occ))
for i in range(occ-1):
    eneroc[i,i+1]=np.sqrt(i+1)
    eneroc[i+1,i]=eneroc[i,i+1]
#Construct number operator as TT
for k in range(qmodes):
    if k==0:
        tmp=tt.kron(tt.matrix(numoc)*freq[k],tt.eye(occ,qmodes-1))
    elif 0<k<qmodes-1:
        tmp=tt.kron(tt.eye(occ,k-1),tt.matrix(numoc)*freq[k])
        tmp=tt.kron(tmp,tt.eye(occ,qmodes-k))
    else:
        tmp=tt.kron(tt.eye(occ,k),tt.matrix(numoc)*freq[k])
    tt_numoc=tt_numoc+tmp
    tt_numoc=tt_numoc.round(eps)
tt_systemnumoc=tt.kron(tt_Ie,tt_numoc)
tt_systemnumoc=tt.kron(tt_systemnumoc,tt.eye(occ,qmodes))
#create a duplicate of number operator for the ficticious system
tt_tildenumoc=tt.kron(tt_Ie,tt.eye(occ,qmodes))
tt_tildenumoc=tt.kron(tt_tildenumoc,tt_numoc)
#initialize displacement operator
tt_energy=tt.eye(occ,qmodes)*0.
tt_tilenergy=tt.eye(occ,qmodes)*0.
for k in range(qmodes):
    if k==0:
#coshtheta takes account for energy flow from real to ficticious system
#thus takes account for temperature effect
        tmp=tt.kron(tt.matrix(eneroc)*gk[k]*coshthetak[k],tt.eye(occ,qmodes-1))
    elif 0<k<qmodes-1:
        tmp=tt.kron(tt.eye(occ,k-1),tt.matrix(eneroc)*gk[k]*coshthetak[k])
        tmp=tt.kron(tmp,tt.eye(occ,qmodes-k))
    else:
        tmp=tt.kron(tt.eye(occ,k),tt.matrix(eneroc)*gk[k]*coshthetak[k])
    tt_energy=tt_energy+tmp
    tt_energy=tt_energy.round(eps)
tt_systemenergy=tt.kron(tt.matrix(pz),tt_energy)
tt_systemenergy=tt.kron(tt_systemenergy,tt.eye(occ,qmodes))
for k in range(qmodes):
    if k==0:
        tmp=tt.kron(tt.matrix(eneroc)*gk[k]*sinhthetak[k],tt.eye(occ,qmodes-1))
    elif 0<k<qmodes-1:
        tmp=tt.kron(tt.eye(occ,k-1),tt.matrix(eneroc)*gk[k]*sinhthetak[k])
        tmp=tt.kron(tmp,tt.eye(occ,qmodes-k))
    else:
        tmp=tt.kron(tt.eye(occ,k),tt.matrix(eneroc)*gk[k]*sinhthetak[k])
    tt_tilenergy=tt_tilenergy+tmp
    tt_tilenergy=tt_tilenergy.round(eps)
tt_tildeenergy=tt.kron(tt.matrix(pz),tt.eye(occ,qmodes))
tt_tildeenergy=tt.kron(tt_tildeenergy,tt_tilenergy)
#Note that ficticious Harmonic oscillators carry negative sign
H=tt_He+tt_systemnumoc-tt_tildenumoc+tt_systemenergy+tt_tildeenergy
H=H.round(eps)
#Construct propagation operator, d/dt psi(t0)=A psi(t0) 
A=-EYE*H
y0=tt_psi0 #Initialize wavefunction
#Heaviside functions, for selecting electronic states from overall wavefunction
tt_heavu=tt.kron(tt_su,tt.ones(occ,dim*2))
tt_heavd=tt.kron(tt_sd,tt.ones(occ,dim*2))
#Propagation time step and range
t=np.arange(0,nsc*tau,tau)
t=t*au2ps
#Add noise, for higher rank KSL propagation
radd = 19
#radd = np.array([1,9,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,9,1]) 
radd=np.array([1,9])
radd=np.append(radd,np.repeat(19,qmodes*2-3))
radd=np.append(radd,np.array([9,1]))
#if ( radd > 0 ):
tt_rand=tt.rand(occ,qmodes*2,radd)
#tt_rand=tt_rand*10**(-140)
#for i in range(radd-1):
#    tt_rand=tt_rand+tt.rand(occ,qmodes*2,1)
wxw=tt_rand.to_list(tt_rand)
for i in range(qmodes*2):
    wxw[i]=wxw[i]*0.005
tt_rand=tt_rand.from_list(wxw)
#tt_rand=tt_rand*tt_rand.norm()**(-1) #Renormalize noise
tt_rand=tt.kron(tt.ones(2,1),tt_rand)
#tt_rand=tt.kron(tmp,tt_rand)
y0 = y0+tt_rand*1e-10 #Ensure noise is small
ul=np.array([[1,0],[0,0]])
ur=np.array([[0,1],[0,0]])
tt_ul=tt.matrix(ul)
tt_ur=tt.matrix(ur)
for i in range(qmodes*2):
    tt_ul=tt.kron(tt_ul,tt.eye(occ,1))
    tt_ur=tt.kron(tt_ur,tt.eye(occ,1))
#Initalize population arrays
psu=np.zeros((nsc))
psd=np.zeros((nsc))
coh12=np.zeros((nsc),dtype=complex)
coh21=np.zeros((nsc),dtype=complex)
ptot=np.zeros((nsc))
expec_H=np.zeros((nsc))
expec_Hu=np.zeros((nsc))
#Propagation loop
for ii in range(nsc):
    y0=tt.ksl.ksl(A,y0,tau)
    print(t[ii])
    psu[ii]=np.abs(tt.dot(tt_heavu*y0,tt_heavu*y0))
    psd[ii]=np.abs(tt.dot(tt_heavd*y0,tt_heavd*y0))
    coh12[ii]=tt.dot(tt.matvec(tt_ul,y0),tt.matvec(tt_ur,y0))
    coh21[ii]=tt.dot(tt.matvec(tt_ur,y0),tt.matvec(tt_ul,y0))
    expec_H=tt.dot(y0,tt.matvec(H,y0))
    expec_Hu=tt.dot(tt_heavu*y0,tt.matvec(H,tt_heavu*y0))


#Plot population difference    
#plt.figure(dpi=600)
#plt.xlim(0.,1.)    
#plt.ylim(-1.,1.)             
#plt.xlabel('time(ps)')
#plt.ylabel('Populations')
#plt.plot(t,psu-psd,label='Model 1, sigma_z')
#plt.legend()                     
    
    
np.savetxt('psu.npy',psu)
np.savetxt('psd.npy',psd)
np.savetxt('t.npy',t)
np.savetxt('coh12.npy',coh12)    
np.savetxt('coh21.npy',coh21)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    