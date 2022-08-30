import math
import scipy.special

# Script to generate bath mode parameters for ohmic spectral density with the
#      spin-boson model
# Parameters based on Kelly, Brackbill, Markland JCP 142 (2015).
# See Fig. 2
#
# Kelly et al. did not specify which protocol they followed
# to generate the bath modes. Could be equidistant.
# But for Ohmic Spectral Density, the protocol described in 
# Kim, Kapral, JCP 123 (2005) seems to be standard.
# However this requires an omega_max not provided in the Kelly paper.
# (Equidistant spacing would require an upper limit too.)
#
# Created by Alex Schubert
#
# Most recent modification 8/29/22 by Ellen Mulvihill

#Required for the script
nu        = 60          # number of modes
xi        = 0.4         # dimensionless Kondo parameter (friction)
omega_c   = 2.          # cut-off frequency
alpha     = 0.95        # threshold for omega_max, 95% corresponds to omega_max=10 for omega_c=2 Delta.
tol       = 1e-6        # tolerance

print('Generating bath modes from Ohmic spectral density using Kim/Kapral(2005) discretization.')
print('Units in terms of electronic coupling energy.')
print('Number of modes =', nu)
print('Kondo parameter =', xi)
print('cut-off freq.   =', omega_c)
z = (alpha-1.)/math.exp(1)
omega_temp = (-scipy.special.lambertw(z,-1,tol)-1.)*omega_c
omega_max = math.ceil(omega_temp.real)
print('Omega_max set to:', omega_max)
print('Required omega_N to reach', alpha, 'Quantile: ',  omega_temp)


outputfile1 = open('specdens_parameters/omega_nm_xi%3.1fwc%3.1f_wmax%1d_dofn%2d.txt' %(xi, omega_c, omega_max, nu), 'w')
outputfile4 = open('specdens_parameters/spectral_density_nm_xi%3.1fwc%3.1f_wmax%1d_dofn%2d.txt' %(xi, omega_c, omega_max, nu), 'w')
outputfile5 = open('specdens_parameters/acc_sum_nm_xi%3.1fwc%3.1f_wmax%1d_dofn%2d.txt' %(xi, omega_c, omega_max, nu), 'w')
outputfile6 = open('specdens_parameters/C_nm_xi%3.1fwc%3.1f_wmax%1d_dofn%2d.txt' %(xi, omega_c, omega_max, nu), 'w')

Delta_omega = (omega_c/nu)*(1-math.exp(-omega_max/omega_c))

omega = []
c = []

#Non-equilibriums shift is set to zero for this case.
#Req is defined such that the Hamiltonian is written as 0.5*omega^2(R*1 - 0.5 Req sigma_z)^2
#Then, c = 0.5 omega^2 Req, and the Hamiltonian coupling term is -c*R*sigma_z

for j in range(0,nu):
    omega.append(float(-omega_c*math.log(1.-(j+1)*Delta_omega/omega_c)))
    c.append(float(math.sqrt(xi*Delta_omega)*omega[j]))
    outputfile1.write("%s\n"%omega[j])
    outputfile6.write("%s\n"%c[j])
    outputfile4.write("%s\t%s\n"%(omega[j], ((math.pi/2.)*(c[j]**2)/(omega[j]))))
  
IntJ = []
accsum = 0
l = 0

for k in range(0,20001):
    x = 0+k*(omega_max/15000.)
    if (l < nu):
        if (x > omega[l]):
            accsum += (math.pi/2.)*xi*Delta_omega*omega[l] #accumulated sum
            l += 1
    if (l == nu):
        print('All modes added to accumulated sum.')
        l += 1
    IntJ.append(float(accsum))
    # cumulative distribution function for comparison of the discrete set of
    # modes to the continues distribution plot 1:2 and 1:3 together.
    outputfile5.write("%s\t%s\t%s\n"%(x, IntJ[k], (math.pi/2.)*omega_c*omega_c*xi*((-x/omega_c-1)*math.exp(-x/omega_c)+1)))
  
print('Script finished.')
