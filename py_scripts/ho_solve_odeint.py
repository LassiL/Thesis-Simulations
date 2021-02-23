#%%
import numpy as np
import matplotlib.pyplot as plt
from math import pi
#LSODA as implemented in odeint seems to work the best
from scipy.integrate import odeint

#%% FUNDAMENTAL CONSTANTS
a_latt = 6.479E-10
echarge = 1.60217662E-19 #coulomb
eps_0 = 8.8541878128E-12 #[F*m^{-1} = C/(Nm^2)]
bohr = 0.529E-10 #meter
a_0 = bohr
amu = 1.660E-27 #kg

#%% ELECTRIC FIELD PARAMETERS - VARY THESE
# Need to consider these inside the material though...
# Impedance - skin depth - absorption - pulse broadening

E_0tot = 2E8 # 2E8 gives around 20% modulation in 2,-2,2 reflection if the freq.
#E_0tot = 2.5E10 # = 250 MV/cm should give equal nlo and lo 

#Polarization in (2,-2,2) direction
E_0x = (1/np.sqrt(3)) * E_0tot
E_0y = -(1/np.sqrt(3)) * E_0tot 
E_0z = (1/np.sqrt(3)) * E_0tot

# FWHM of the Gaussian laser pulse: E(t) = (\vec{E_0})*\exp((t/\tau)^2)
#FWHM = 300E-15 #300 fs
FWHM = 300E-15 
tau = FWHM/np.sqrt(2*np.log(2))
alpha = 1/tau**2

#%% Simulated material constants 

lifetime = 5E-12 # UNITS [s] i.e., 10 ps
#According to mahan eq. 7.47 the lifetime is approx 1 ps at 300 K
#This damping is half of the rate from fermi's golden 1/s units
Gamma = (1/lifetime)/2 #Damping [1/s], note Gamma in my H.O eqns multiplied by 2
#Gamma = 1061136396569.9978/2 #From mahan.py
omega_0 = 5.76E12 * 2*pi # Ph_2.out

#Primitive cell volume m^3
Omega_0 = 431 * bohr**3 

#Indium 
m_In = 114.82 * amu #mass [kg]
#Z_In = 1.94 * echarge # Born effective [Coulomb]
Z_In = 2.05 * echarge # Ph_2.out
d_In = -0.43 * (4*np.pi*eps_0)/a_0 #Raman tensor 


#%% Define time interval globally
# This needs to be imported to trxd.py script together with the time steps.
# make it a function here and call in txrd
ts = np.linspace(0, 10E-12, 10000)

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

plt.figure(1)
#plt.title("E_0tot = 2E8 V/m - x-direction")
plt.xlabel("time (ps)")
plt.ylabel("$\Delta$x (pm)")
plt.rcParams.update({'font.size': 20})
plt.grid()
plt.tight_layout()
plt.plot(ts, Uxs[:,0])

plt.xlim(0.5E-12, max(ts))
plt.xticks([0.2E-11, 0.4E-11, 0.6E-11, 0.8E-11, max(ts)], (2, 4, 6, 8, 10)) #[pm]
#For 20 ps range
#plt.xticks([0.5E-11, 1.0E-11, 1.5E-11, max(ts)])
#plt.xticks(np.arange(0.5E-12,2E-11,0.5E-11)) 

plt.ylim(-1.6E-13, 1.6E-13)
#plt.yticks([-1.0E-13, -0E-13, 1.0E-13], (-0.1, 0, 0.1)) #[pm]
plt.yticks([-1.5E-13, -1.0E-13, -0.5E-13, -0E-13, 0.5E-13, 1.0E-13, 1.5E-13], (-0.15, '0.10', 0.05, 0, 0.05, '0.10', 0.15))

#For 20 ps plot
#plt.yticks(np.arange(-1.5E-13, 1.55E-13, 0.5E-13))

#plt.tight_layout()
plt.savefig('figures/harm_osc/2E8_x.png', format='png', dpi=300)
plt.show()

#%% Figure for the proposal

plt.figure(2)

plt.xlabel("time (ps)")
#What do i label this
plt.ylabel("$\Delta \overrightarrow{r}_{In}$ (pm)")
plt.rcParams.update({'font.size': 20})
plt.grid()
#Multiply by sqrt(3) to make it the horizontal direction
plt.tight_layout()
plt.plot(ts, np.sqrt(3)*Uxs[:,0])

plt.xlim(0.565E-12, max(ts))
plt.xticks([0.2E-11, 0.4E-11, 0.6E-11, 0.8E-11, max(ts)], (2, 4, 6, 8, 10)) #[pm]
#For 20 ps range
#plt.xticks([0.5E-11, 1.0E-11, 1.5E-11, max(ts)])
#plt.xticks(np.arange(0.5E-12,2E-11,0.5E-11)) 

plt.ylim(-3.2E-13, 3.2E-13)
plt.yticks([-3.0E-13, -2.0E-13, -1.0E-13, 0, 1E-13, 2E-13, 3E-13], 
(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))
#plt.yticks([-1.0E-13, -0E-13, 1.0E-13], (-0.1, 0, 0.1)) #[pm]
#plt.yticks([-2.5E-13, -2.0E-13, -1.5E-13, -1.0E-13, -0.5E-13, -0E-13, 0.5E-13, 1.0E-13, 1.5E-13, 2.0E-13, 2.5E-13], (-0.25, '-0.20', -0.15, '-0.10', -0.05, 0, 0.05, '0.10', 0.15, '0.20', 0.25))

#For 20 ps plot
#plt.yticks(np.arange(-1.5E-13, 1.55E-13, 0.5E-13))

plt.tight_layout()
plt.savefig('figures/harm_osc/2E8_x.png', format='png', dpi=300)
plt.show()

#%% For Proposal PLOT delta a_bond ! essentially two times the above

plt.figure(3)

plt.xlabel("time (ps)")
#What do i label this
plt.ylabel("$\Delta a_{bond}$ (pm)")
plt.rcParams.update({'font.size': 20})
plt.grid()
#Multiply by sqrt(3) to make it the horizontal direction

plt.plot(ts, 2*np.sqrt(3)*Uxs[:,0])

plt.xlim(0.565E-12, max(ts))
plt.xticks([0.2E-11, 0.4E-11, 0.6E-11, 0.8E-11, max(ts)], (2, 4, 6, 8, 10)) 

plt.ylim(-0.6E-12, 0.6E-12)
plt.yticks([-5E-13, -2.5E-13, 0, 2.5E-13, 5E-13], ('-0.50', -0.25, 0, 0.25, '0.50'))
# plt.yticks([-3.0E-13, -2.0E-13, -1.0E-13, 0, 1E-13, 2E-13, 3E-13], 
# (-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))



plt.tight_layout()
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

plt.figure(4)
#plt.title("E_0tot = 2E8 V/m - y-direction")
plt.xlabel("time (s)")
plt.ylabel("dy (m)")
plt.plot(ts, Uys[:,0])
plt.savefig('figures/harm_osc/2E8_y.png', format='png', dpi=300)
#plt.tight_layout()
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

plt.figure(5)
#plt.title("E_0tot = 2E8 V/m - z-direction")
plt.xlabel("time (s)")
plt.ylabel("dz (m)")
plt.plot(ts, Uzs[:,0])
plt.tight_layout()
plt.savefig('figures/harm_osc/2E8_z.png', format='png', dpi=300)
plt.show()

# #%% This bond length change is kinda useless
# Plotting the change in the bond length da_bond/a_bond
# Oscillation along the bond length
# Assuming that In and Sb move equal amounts but in the opposite
# Directions
a_bond = 2.806E-10
#Factor of two from In Sb moving against each other
da_bond = 2*np.sqrt(np.square(Uxs[:,0])+np.square(Uys[:,0])+np.square(Uzs[:,0]))
#This does not work.......
daa_bond = np.linspace(-da_bond[0],da_bond[0], 10000)
plt.figure(6)
plt.plot(ts, da_bond/a_bond*100)
#plt.xlim(0.55E-12,max(ts))

plt.xlim(0.565E-12, max(ts))
plt.xticks([0.2E-11, 0.4E-11, 0.6E-11, 0.8E-11, max(ts)], (2, 4, 6, 8, 10)) 
plt.ylim(0, 0.2)
plt.ylabel("$\Delta a_{bond}/a_{bond}$ (%)")
plt.grid()
plt.xlabel("time (ps)")
plt.tight_layout()
plt.savefig('figures/harm_osc/da_bond.png', format='png', dpi=300)
plt.show()


# %%
