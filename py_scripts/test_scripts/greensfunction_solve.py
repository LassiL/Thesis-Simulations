import numpy as np
import matplotlib.pyplot as plt
"""THIS SCRIPT IS FOR COMPARING ODEINT (LSODA) AND INTEGRATION OF THE
GREEN'S FUNCTION SOLUTION"""
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
omega_0 = 5.85E12 #Frequency Hz Should this be angular or f?
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

#%%
"""ODEINT For reference"""
from scipy.integrate import odeint
#x-direction
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

Ux0 = [0,0]
ts = np.linspace(0,50E-12,50000)
Uxs = odeint(dUx_dt, Ux0, ts)

plt.figure()
plt.title("E_0tot = 1E7 V/m - ODEINT (def. LSODA)")
plt.xlabel("time (s)")
plt.ylabel("dx (m)")
plt.plot(ts, Uxs[:,0])
#plt.savefig('1E7_x.eps', format='eps', dpi=300)
plt.show()

# #%%
# "EXAMPLE FROM: https://www.southampton.ac.uk/~fangohr/teaching/python/book/html/16-scipy.html"

# from math import cos, exp, pi
# from scipy.integrate import quad

# # function we want to integrate
# def f(x):
#     return exp(cos(-2 * x * pi)) + 3.2

# # call quad to integrate f from -2 to 2
# res, err = quad(f, -2, 2)

# print("The numerical result is {:f} (+-{:g})"
#     .format(res, err))

#%%
"""GREEN'S FUNCTION - this works and save for reference"""

def x(t_pr): #What is this function of? t_pr is the integration variable
    #But it is actually a function of t
    #Define t outside the function
    t = 2E-12 #This is both a constant in the integrand and in the upper limit
    A = exp(-Gamma*t)/(m_In*sqrt(omega_0**2 - Gamma**2))
    c = sqrt(omega_0**2 - Gamma**2)
    a = Z_In * E_0x
    b = 2 * Omega_0 * d_In * E_0y * E_0z
    #print(A,c,a,b)
    #probably dont need np.squares here
    F_x = a*np.exp(-alpha*t_pr**2)+b*np.exp(-2*alpha*t_pr**2)
    #F_x = (a*exp(-alpha*np.square(t_pr)+b*exp(-2*alpha*np.square(t_pr))))
    #Note the t in return is also the limit value I am integrating
    #Need to make a linspace of it
    #And t_pr is what I am integrating respect to
    return A*sin(c*(t-t_pr))*exp(Gamma*t_pr)*F_x


#Prints the result and the error
sol_ref = quad(x, -2E-12, 2E-12)
print("The numerical reference result is {}".format(sol_ref))

#%%
"Green's function - need to make the t to be evaluated over an interval"
#How to define the t both in function and in the limit as an interval

def x(t_pr): #What is this function of? t_pr is the integration variable
    #But it is actually a function of t
    #Define t outside the function
    #t = 2.2E-12 #This is both a constant in the integrand and in the upper limit
    #A = exp(-Gamma*t)/(m_In*sqrt(omega_0**2 - Gamma**2))
    c = sqrt(omega_0**2 - Gamma**2)
    a = Z_In * E_0x
    b = 2 * Omega_0 * d_In * E_0y * E_0z
    #print(A,c,a,b)
    #probably dont need np.squares here
    F_x = a*np.exp(-alpha*t_pr**2)+b*np.exp(-2*alpha*t_pr**2)
    #F_x = (a*exp(-alpha*np.square(t_pr)+b*exp(-2*alpha*np.square(t_pr))))
    #Note the t in return is also the limit value I am integrating
    #Need to make a linspace of it
    #And t_pr is what I am integrating respect to
    return sin(c*(t-t_pr))*exp(Gamma*t_pr)*F_x

#t = np.linspace(2E-12, 2.5E-12, 5)
t = 2E-12
#Need to make A an empty list and multiply the quad with it
#A = []
A = exp(-Gamma*t)/(m_In*sqrt(omega_0**2 - Gamma**2))

# for i in range(len(t)):
#     print(A[i]*quad(x, -2E-2, t[i]))

sol = quad(x, -2E-12, t)

print("The numerical result is {}".format(A*sol[0]))
# %%
