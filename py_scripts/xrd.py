#%% Loading packages
import numpy as np 
import matplotlib.pyplot as plt 
#from math import sqrt, sin, asin, pi 
import math

#%% Constants 
a = 6.479E-10 # Lattice constant of InSb at 300 K (Landolt-Börnstein) [m]
#a = 6.479 # [Å]

#Planck constant (unreduced) times speed of light [eV*m]
h_Planck = 4.135667696E-15 # Planck constant (unreduced) [eV*s]
c = 3.00E8 # Speed of light [m/s]
#c = 3.00E18 # Å/s
hc = h_Planck * c # [eV*m]

#%% Spacing of different crystal planes
def d_hkl(h,k,l):
    return a/sqrt(h**2 + k**2 + l**2)

print(d_hkl(2,2,2))

# %% Bragg's law respect to theta_hkl
def theta_hkl(h,k,l):
    n = 1 # Diffraction order
    return 180/pi * asin(n*hc/(E*2*d_hkl(h,k,l)))

#%% Make plots of the Bragg angles for an energy interval
# Best would probably be to create a slit of h,k,l values...
# Or include this in the above function

Energy = np.linspace(3.56, 10, 1000)*1000 # Available energy range at FemtoMAX

# Interested reflections
# "Forbidden" 
theta_200 = []; theta_222 = []; theta_600 = []; 
# Comparison to the highest intensity
theta_111 = []

for E in Energy:
    theta_200.append(theta_hkl(2,0,0))
    theta_222.append(theta_hkl(2,2,2))
    theta_111.append(theta_hkl(1,1,1))

#600 error with the math domain - goes over 90 deg below 5.8 keV
Energy600 = np.linspace(5.8, 10, 1000)*1000
for E in Energy600:
    theta_600.append(theta_hkl(6,0,0))

# Data tip and axis i.e., fig, ax plt.subplots(1)
plt.figure(1)
plt.plot(Energy, theta_200, label='(200)')
plt.plot(Energy, theta_222, label='(222)')
#plt.plot(Energy, theta_2_22, label='(2-22)')
plt.plot(Energy, theta_111, label='(111)')
plt.plot(Energy, theta_600, label='(600)')
#plt.title('title')
plt.xlabel('E (eV)')
plt.ylabel(r'$\theta$ (deg)')
plt.legend()
plt.grid()
plt.savefig('figures/theta(E,h,k,l).png', format='png', dpi=300)
plt.show()

# %%

# %%
