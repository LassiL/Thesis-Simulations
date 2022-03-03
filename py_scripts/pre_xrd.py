#%%
"""This script calculates the scattering vectors of the different nearly forbidden reflections 
and their atomic form factors to be used in the txrd.py script"""

#%% Loading packages
import numpy as np 
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size 
from math import sqrt, pi 

#%% Define constants 
a = 6.479E-10 # Lattice constant of InSb at 300 K (Landolt-Börnstein) [m]

h_Planck = 4.135667696E-15 # Planck constant (unreduced) [eV*s]
c = 3.00E8 # Speed of light [m/s]
hc = h_Planck * c  

#%% Spacing of different crystal planes
def d_hkl(h,k,l):
    return a/sqrt(h**2 + k**2 + l**2)

#%% Scattering vectors for the based on the above bragg angles  
def q_hkl(h,k,l):
    return 2*pi/(d_hkl(h,k,l)*1E10) #In Ångströms

#%% Atomic form factor for the structure factor as parametrized in
# http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php and reference therein (Journal of crystallography)

#Indium atom
def f_In(q):
    # Constants for the Gaussian parametrization of the atomic form factor
    a_In = np.array([19.1624, 18.5596, 4.2948, 2.0396])
    b_In = np.array([0.5476, 6.3776, 25.8499, 92.8029])
    c_In = 4.9391
    return np.sum(a_In*np.exp(-1*b_In*(q/(4*pi))**2)) + c_In

#Antimony atom 
def f_Sb(q):
    # Constants for the Gaussian parametrization of the atomic form factor
    a_Sb = np.array([19.6418, 19.0455, 5.0371, 2.6827])
    b_Sb = np.array([5.3034, 0.4607, 27.9074, 75.2825])
    c_Sb = 4.5909
    return np.sum(a_Sb*np.exp(-1*b_Sb*(q/(4*pi))**2)) + c_Sb

#%% Create arrays for the form factors for plotting
f_In1 = []; f_Sb1 = []

interval = np.linspace(0,8,10000) #define the scattering vector interval
for q in interval:
    f_In1.append(f_In(q))
    f_Sb1.append(f_Sb(q))

f_In_arr = np.array(f_In1)
f_Sb_arr = np.array(f_Sb1)

#%% Plot the form factors
plt.figure(555)
#plot the form factors for the chosen q-vector interval 
plt.plot(interval, f_In_arr, label="$f_{In}$")
plt.plot(interval, f_Sb_arr, label="$f_{Sb}$")
#plot the scattering vectors of the different reflections as vertical lines
x_ticks = [0, 4, 8, q_hkl(2,0,0), q_hkl(2,2,2), q_hkl(6,0,0)] 
labels = [0, 4, 8, r'$q_{200}$', r'$q_{222}$', r'$q_{600}$']
plt.xticks(x_ticks, labels, rotation ='horizontal')

#Some formatting of the plot
plt.xlabel(r'$q$ (Å$^{-1}$)', size='14')
#r'$\theta$ (deg)'
plt.ylabel(r'$f(q)$', size='14')
plt.legend()
plt.xlim(0, 8)
plt.ylim(15, 55)
#plt.ylim(0, max(f_In_arr))
plt.tight_layout()
plt.grid()
#plt.savefig('figures/pre_xrd/form_factors.png', format='png', dpi=300)
plt.show()

# %%
