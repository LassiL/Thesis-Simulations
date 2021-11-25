#%%
"""This scripts solves the driven damped harmonic oscillators for the user defined THz pulse parameters 
and anharmonic lifetime given the DFPT simulated material parameters"""

#%% Lodaing packages
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from numpy.lib.type_check import real

# LSODA algorithm for solving the harmonic oscillators 
# as implemented in odeint seems to work the best
from scipy.integrate import odeint

"""Definition of the THz laser pulse"""
#%% THz field parameters for the driving force of the harmonic oscillator 
# VARY THESE PARAMETERS

#####DEFINE THE POLARIZATION AND THE FIELD STRENGTHS#########
#Total field strength
E_0tot = 1.00E8 # in units of V/m 

#Polarization in (2,-2,2) direction
#E_0x = (1/np.sqrt(3)) * E_0tot
#E_0y = -(1/np.sqrt(3)) * E_0tot 
#E_0z = (1/np.sqrt(3)) * E_0tot

##Polarization only in x-direction
E_0x = E_0tot
E_0y = 0
E_0z = 0

######DEFINE THE PULSE PARAMETERS########
# FWHM of the Gaussian laser pulse: E(t) = (\vec{E_0})*\exp((t/\tau)^2)
FWHM = 200E-15 # seconds
#Define for convienience
tau = FWHM/np.sqrt(2*np.log(2))
alpha = 1/tau**2

#Carrier frequency
omega_THz = 2.8*np.pi*3E12

#Center the pulse at 1 ps for easier numerical treatment
t_0 = 1E-12 #center of the gaussian pulse

#Visualize the pulse
plt.figure(222)
ttt = np.linspace(0.5E-12, 1.5E-12, 20000)

plt.plot(ttt, np.real(E_0tot*np.exp(-alpha*np.square(ttt-t_0))*np.exp(1j*omega_THz*(ttt-t_0))),label='E-field')
plt.plot(ttt, np.real(E_0tot*np.exp(-alpha*np.square(ttt-t_0))), label='Envelope' ) 
plt.legend()
plt.tight_layout()
plt.xlim(0.5E-12,1.5E-12)
plt.xlabel("t (s)")
plt.ylabel(r"$E(t)$ (V/m)")
plt.show()

#%% Defining fundamental constants
a_latt = 6.479E-10
echarge = 1.60217662E-19 #coulomb
eps_0 = 8.8541878128E-12 #[F*m^{-1} = C/(Nm^2)]
bohr = 0.529E-10 #meter
a_0 = bohr
amu = 1.660E-27 #kg

#%%  Anharmonc lifetime and the DFPT simulated InSb material parameters 

## Anharmonic lifetime defined in the harmonic oscillators as
# 2gamma = 1/tau  
# Ferry: tau = 7.26E-12 [seconds]
# Lockwood dielectric 2gamma = 3.37 cm-1 
# Palik dielectric 2gamma = 2.86 cm-1
gamma = (2.085 * 2.998E10) # 8 ps lifetime

#LO mode angular frequency 
omega_0 = 5.77E12 * 2*pi # [2\pi THz] 

#Primitive cell volume m^3
Omega_0 = 1/4 * a_latt**3 # experimental 

#Indium 
m_In = 114.82 * amu #mass [kg]
Z_In = 2.03 * echarge # simulated Born effective charge 
d_In = -0.559845 * (4*np.pi*eps_0)/a_0 # simulated Raman tensor 

#%%
"""Solving the harmonic oscillator equations in x y and z directions for the above defined THz pulse that 
will be then imported to the txrd script to plot the time resolved x-ray diffraction of the chosen reflections 
(we will only solve them for displacements of the Indium atoms as we can scale them later for Sb)"""

#%% Define the time interval to be used for solving the harmonic oscillators 
#Note: THE VALUE USED HERE MUST MATCH THE VALUE USED IN THE txrd.py SCRIPT 
ts = np.linspace(0, 11E-12, 40000)

#%% Harmonic oscllator for the x-direction
def dUx_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    a = Z_In * E_0x
    b = 2 * Omega_0 * d_In * E_0y * E_0z
    #F_x = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0))) #only LO
    #F_x = b*np.exp(-2*alpha*np.square(t-t_0))*np.real(np.exp(2*1j*omega_THz*(t-t_0))) #only NLO
    F_x = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0)))+b*np.exp(-2*alpha*np.square((t-t_0)))*np.real(np.exp(2*1j*omega_THz*(t-t_0)))
    return [U[1], F_x/m_In - 2*gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0] #Initial conditions
Uxs = odeint(dUx_dt, U0, ts) #Solve the harmonic oscillators in x direction
#Save the displacements into a text file to be later imported to the txrd script
np.savetxt("data/dirx.txt", Uxs[:,0])

#Plot the displacements
plt.figure(1)
plt.xlabel("time (s)")
plt.ylabel("$u_{In,x}$ (m)")
plt.rcParams.update({'font.size': 20})
plt.grid()
plt.tight_layout()
plt.plot(ts, Uxs[:,0])
plt.xlim(0, max(ts))
#plt.ylim(-2E-13, 2E-13)
#plt.savefig('figures/harm_osc/1E8_x.png', format='png', dpi=300)
plt.show()

#%% Harmonic oscillator in the y-direction
def dUy_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    a = Z_In * E_0y
    b = 2 * Omega_0 * d_In * E_0x * E_0z
    #F_y = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0))) #LO
    #F_y = b*np.exp(-2*alpha*np.square(t-t_0))*np.real(np.exp(2*1j*omega_THz*(t-t_0))) #NLO
    F_y = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0)))+b*np.exp(-2*alpha*np.square((t-t_0)))*np.real(np.exp(2*1j*omega_THz*(t-t_0)))
    return [U[1], F_y/m_In - 2*gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0]
Uys = odeint(dUy_dt, U0, ts)
np.savetxt("data/diry.txt", Uys[:,0])

plt.figure(2)
plt.xlabel("time (s)")
plt.ylabel("$u_{In,y}$ (m)")
plt.rcParams.update({'font.size': 20})
plt.grid()
plt.tight_layout()
plt.xlim(0, max(ts))
#plt.ylim(-2E-13, 2E-13)
plt.plot(ts, Uys[:,0])
#plt.savefig('figures/harm_osc/1E8_y.png', format='png', dpi=300)
plt.show()

#%% Harmonic oscillator in the z-direction

def dUz_dt(U, t):
    # Vector U with x=U[0] and v_x=U[1] and the function returns
    # [x' and v_x']
    a = Z_In * E_0z
    b = 2 * Omega_0 * d_In * E_0x * E_0y
    #F_z = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0))) #LO
    #F_z = b*np.exp(-2*alpha*np.square(t-t_0))*np.real(np.exp(2*1j*omega_THz*(t-t_0))) #NLO
    F_z = a*np.exp(-alpha*np.square(t-t_0))*np.real(np.exp(1j*omega_THz*(t-t_0)))+b*np.exp(-2*alpha*np.square((t-t_0)))*np.real(np.exp(2*1j*omega_THz*(t-t_0)))
    return [U[1], F_z/m_In - 2*gamma*U[1] - omega_0**2*U[0]]

U0 = [0,0]
Uzs = odeint(dUz_dt, U0, ts)
np.savetxt("data/dirz.txt", Uzs[:,0])

plt.figure(3)
plt.xlim(0, max(ts))
#plt.ylim(-2E-13, 2E-13)
plt.xlabel("time (s)")
plt.ylabel("$u_{In,z}$ (m)")
plt.rcParams.update({'font.size': 20})
plt.grid()
plt.tight_layout()
plt.plot(ts, Uzs[:,0])
plt.tight_layout()
#plt.savefig('figures/harm_osc/1E8_z.png', format='png', dpi=300)
plt.show()
# %%
