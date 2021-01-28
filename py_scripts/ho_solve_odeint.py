#%%
import numpy as np
import matplotlib.pyplot as plt
from math import pi
#LSODA as implemented in odeint seems to work the best
from scipy.integrate import odeint

#%% Constants
# Fundamental 
a_latt = 6.479E-10

echarge = 1.60217662E-19 #coulomb
eps_0 = 8.8541878128E-12 #[F*m^{-1} = C/(Nm^2)]
bohr = 0.529E-10 #meter
a_0 = bohr
amu = 1.660E-27 #kg

# Simulated
damping = 10E-12 # this damping works with the 5.85E12 frequency
#damping = 40E-12*2*np.pi
Gamma = 1/damping #Damping seconds - Note Gamma = 1/(damping_rate) = 1/1E-12
#Gamma = 0.03E12
#omega_0 = 5.85E12 #Frequency Hz THIS IS ACTUALLY THE FREQUENCY f not omega ...
# THE ACTUAL ONE IS THIS
omega_0 = 5.85E12 * 2*pi # 36.76E12 rad/s 
# Solution becomes extremely unstable if i multiply the frequency by 2pi ...
# Need to try different solvers
"""In the diamond papers they use 40 Thz and 0.3 ps-1 damping (oscillations 
decay in less than 10 ps)"""
"""BUT I THINK IT IS ACTUALLY THE frequency (not omega) as diamond optical mode is 
1300cm-1 -> f = 1300 * 3E10 = 39 THz"""

#Define time interval globally
ts = np.linspace(0, 20E-12, 10000)

Omega_0 = 431 * bohr**3 #Primitive cell volume m^3

#Indium 
m_In = 114.82 * amu #mass [kg]
Z_In = 1.94 * echarge #Born effective [Coulomb]
d_In = -0.43 * (4*np.pi*eps_0)/a_0 #Raman tensor [unit?]

#%%
# Electric fields 
# Need to consider these inside the material though...
# Impedance - skin depth - absorption - pulse broadening

# E_0x = 1/sqrt(3)*E_0tot
#E_0tot = 1E6 # = 10 kV/cm gives 
#50 kV/cm seems to give the right dx = 1.5E-13
#E_0tot = 5E6 # = 50 kV/cm = 0.05 MV/cm gives 
#E_0tot = 1E7 # = 100 kV/cm = 0.1 MV/cm !this gives correct for 5.85 Thz
# Need like 2 MV/cm to get the 
E_0tot = 2E8 # 2E8 gives around 20% modulation in 2,-2,2 reflection if the freq.
# is multiplied by factor of 2pi
#E_0tot = 5E8 # = 5000 kV/cm = 5 MV/cm gives 
#E_0tot = 2.5E10 # = 250 MV/cm should give equal nlo and lo 

#Polarization in (2,-2,2) direction
E_0x = (1/np.sqrt(3)) * E_0tot
E_0y = -(1/np.sqrt(3)) * E_0tot 
E_0z = (1/np.sqrt(3)) * E_0tot

# FWHM of the Gaussian laser pulse: E(t) = (\vec{E_0})*\exp((t/\tau)^2)
FWHM = 50E-15 #300 fs
tau = FWHM/np.sqrt(2*np.log(2))
alpha = 1/tau**2

#%% x-direction
def dUx_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    # Maybe remove the np.square(t)?
    a = Z_In * E_0x
    b = 2 * Omega_0 * d_In * E_0y * E_0z
    #F_x = a*np.exp(-alpha*np.square(t)) #only linear
    #F_x = b*np.exp(-2*alpha*np.square(t)) #only nlo
    F_x = a*np.exp(-alpha*np.square(t))+b*np.exp(-2*alpha*np.square(t))
    return [U[1], F_x/m_In - 2*Gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0]
#ts = np.linspace(0,50E-12,10000)
Uxs = odeint(dUx_dt, U0, ts)

np.savetxt("data/dirx.txt", Uxs[:,0])

#In units of lattice constant
#np.savetxt("data/dirx.txt", Uxs[:,0]/a_latt)

plt.figure()
plt.title("E_0tot = 2E8 V/m - x-direction")
plt.xlabel("time (s)")
plt.ylabel("dx (m)")
plt.plot(ts, Uxs[:,0])
plt.savefig('figures/harm_osc/1E7_x.png', format='png', dpi=300)
plt.show()

#%% y-direction
def dUy_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    # Maybe remove the np.square(t)?
    a = Z_In * E_0y
    b = 2 * Omega_0 * d_In * E_0x * E_0z
    #F_y = a*np.exp(-alpha*np.square(t)) #only linear
    #F_y = b*np.exp(-2*alpha*np.square(t)) #only nlo
    F_y = a*np.exp(-alpha*np.square(t))+b*np.exp(-2*alpha*np.square(t))
    return [U[1], F_y/m_In - 2*Gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0]
#ts = np.linspace(0,50E-12,10000)
Uys = odeint(dUy_dt, U0, ts)

np.savetxt("data/diry.txt", Uys[:,0])

#In units of lattice constant
#np.savetxt("data/diry.txt", Uys[:,0]/a_latt)

plt.figure()
plt.title("E_0tot = 2E8 V/m - y-direction")
plt.xlabel("time (s)")
plt.ylabel("dy (m)")
plt.plot(ts, Uys[:,0])
plt.savefig('figures/harm_osc/1E7_y.png', format='png', dpi=300)
plt.show()

#%% z-direction

def dUz_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    # Maybe remove the np.square(t)?
    a = Z_In * E_0z
    b = 2 * Omega_0 * d_In * E_0x * E_0y
    #F_z = a*np.exp(-alpha*np.square(t)) #only linear
    #F_z = b*np.exp(-2*alpha*np.square(t)) #only nlo
    F_z = a*np.exp(-alpha*np.square(t))+b*np.exp(-2*alpha*np.square(t))
    return [U[1], F_z/m_In - 2*Gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0]
#ts = np.linspace(0,50E-12,10000)
Uzs = odeint(dUz_dt, U0, ts)

np.savetxt("data/dirz.txt", Uzs[:,0])

#In units of lattice constant
#np.savetxt("data/dirz.txt", Uzs[:,0]/a_latt)

plt.figure()
plt.title("E_0tot = 2E8 V/m - z-direction")
plt.xlabel("time (s)")
plt.ylabel("dz (m)")
plt.plot(ts, Uzs[:,0])
plt.savefig('figures/harm_osc/1E7_z.png', format='png', dpi=300)
plt.show()

# #%% This bond length change is kinda useless
# # Plotting the change in the bond length da_bond/a_bond
# # Oscillation along the bond length
# # Assuming that In and Sb move equal amounts but in the opposite
# # Directions
# a_bond = 2.806E-10
# #Factor of two from In Sb moving against each other
# da_bond = 2*np.sqrt(np.square(Uxs[:,0])+np.square(Uys[:,0])+np.square(Uzs[:,0]))
# #This does not work.......
# daa_bond = np.linspace(-da_bond[0],da_bond[0],10000)
# plt.figure()
# plt.plot(ts, da_bond/a_bond*100)
# plt.ylabel("abs(da_bond)/a_bond (%)")
# plt.xlabel("time (s)")
# plt.savefig('figures/harm_osc/da_bond.png', format='png', dpi=300)
# plt.show()

