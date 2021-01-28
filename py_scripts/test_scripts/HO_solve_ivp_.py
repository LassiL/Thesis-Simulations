import numpy as np
import matplotlib.pyplot as plt
"""THIS SCRIPT IS FOR COMPARING METHODS - BUT THE CONCLUSION IS TO USE
THE ODEINT INSTEAD OF SOLVE_IVP"""
# CONSTANTS #
# Fundamental 
echarge = 1.60217662E-19 #coulomb
eps_0 = 8.8541878128E-12 #[F*m^{-1} = C/(Nm^2)]
bohr = 0.529E-10 #meter
a_0 = bohr 
amu = 1.660E-27 #kg

# Simulated
# Note my H.O. eq. has 2*Gamma so the damping is 
damping = 40E-12 #40 ps
Gamma = 1/damping #Damping seconds - Note Gamma = 1/(damping_rate) = 1/1E-12
omega_0 = (2*np.pi)*5.85E12 #Frequency Hz Should this be angular or f?
Omega_0 = 431 * bohr**3 #Primitive cell volume m^3

#Indium 
m_In = 114.82 * amu #mass [kg]
Z_In = 1.94 * echarge #Born effective [Coulomb]
# Indium is positively charged by +2 in the polar crystal.
# Does Sb actually have higher electronegativity?
d_In = -0.43 * (4*np.pi*eps_0)/a_0 #Raman tensor [unit?]

#Electric fields
# Electric fields (Polarization in 2,-2,2 direction)

#E_0tot = 1E6 # = 10 kV/cm gives 
#50 kV/cm seems to give the right dx = 1.5E-13
#E_0tot = 5E6 # = 50 kV/cm = 0.05 MV/cm gives 
E_0tot = 1E7 # = 100 kV/cm = 0.1 MV/cm !USE THIS!!
#E_0tot = 1E8 # = 1000 kV/cm = 1 MV/cm gives 
#E_0tot = 2.5E10 # = 250 MV/cm should give equal nlo and lo 


### NEED TO WRITE A FUNCTION THAT TAKES THE FIELD AS AN INPUT!
#How do I get the 2,-2,2 direction here???????
E_0x = (1/np.sqrt(3)) * E_0tot
E_0y = -(1/np.sqrt(3)) * E_0tot #Mby this is just negative?
E_0z = (1/np.sqrt(3)) * E_0tot

# FWHM of the Gaussian laser pulse
# E(t) = (\vec{E_0})*\exp((t/\tau)^2)
# Should consider how the field actually propagates inside the material...
FWHM = 300E-15 #300 fs
tau = FWHM/np.sqrt(2*np.log(2))
alpha = 1/tau**2

#%% Comment out
# """ODEINT For reference"""
# from scipy.integrate import odeint
# #x-direction
# def dUx_dt(U, t):
#     # Vector U with x=U[0] and v_x=U[1] and the function returns
#     # [x' and v_x']
#     # Maybe remove the np.square(t)?
#     a = Z_In * E_0x
#     b = 2 * Omega_0 * d_In * E_0y * E_0z
#     #F_x = a*np.exp(-alpha*np.square(t)) #only linear
#     #F_x = b*np.exp(-2*alpha*np.square(t)) #only nlo
#     F_x = a*np.exp(-alpha*np.square(t))+b*np.exp(-2*alpha*np.square(t))
#     return [U[1], F_x/m_In - 2*Gamma*U[1] - omega_0**2*U[0]]

# Ux0 = [0,0]
# ts = np.linspace(0,10E-12,50000)
# Uxs = odeint(dUx_dt, Ux0, ts)

# plt.figure()
# plt.title("E_0tot = 1E7 V/m - ODEINT (def. LSODA)")
# plt.xlabel("time (s)")
# plt.ylabel("dx (m)")
# plt.plot(ts, Uxs[:,0])
# #plt.savefig('1E7_x.eps', format='eps', dpi=300)
# #plt.show()

#%% 
"""Solve_IVP"""
# Can input a method, RK45 doesn't work for stiff ODE tho.
# Probably the best is just to use the odeint (LSODA)
from scipy.integrate import solve_ivp
def dUx_dt(t, U):
    a = Z_In * E_0x
    b = 2 * Omega_0 * d_In * E_0y * E_0z
    #F_x = a*np.exp(-alpha*np.square(t)) #only linear
    #F_x = b*np.exp(-2*alpha*np.square(t)) #only nlo
    #Maybe should write the force as a function of E_0x, E_0y, E_0z
    F_x = a*np.exp(-alpha*np.square(t))+b*np.exp(-2*alpha*np.square(t))
    return [U[1], F_x/m_In - 2*Gamma*U[1] - omega_0**2*U[0]]

t_span = [0,10E-12]
#Initial values
x0 = 0
v_x0 = 0
Ux0 = [x0,v_x0]
#Can test different methods, but seems LSODA is the best, which is the 
#Default for odeint... 
#RK23 damps a lot faster than other methods, I think LSODA is the best.
#Radau, BDF or LSODA for stiff problems
solved = solve_ivp(dUx_dt, t_span, Ux0, method='BDF', max_step=1E-16)

plt.figure()
plt.title("E_0tot = 1E7 V/m - RK45")
plt.xlabel("time (s)")
plt.ylabel("dx (m)")
plt.plot(solved.t, solved.y[0]) #Note the solve_ivp always has t and y attributes
plt.show()
# %%
