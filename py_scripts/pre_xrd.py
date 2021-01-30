#%% Loading packages
import numpy as np 
import matplotlib.pyplot as plt 
from math import sqrt, asin, pi 

#import math
#from numpy.core.function_base import linspace 

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

# Should really use a list of h,k,l:
# h = [2,2,1,6]; k = [0,2,1,0]; l = [0, 2, 1, 0]
# print(d_hkl(2,2,2))

# %% Bragg's law respect to theta_hkl
def theta_hkl(h,k,l):
    n = 1 # Diffraction order
    return 180/pi * asin(n*hc/(E*2*d_hkl(h,k,l)))

#%% Make plots of the Bragg angles for an energy interval
# Best would probably be to create a slit of h,k,l values...
# Or include this in the above function

Energy = np.linspace(3.56, 10, 1000)*1000 # Available energy range at FemtoMAX

# Interested reflections:
# "Forbidden" 
theta_200 = []; theta_222 = []; theta_600 = []; 
# Comparison to the highest intensity
theta_111 = []

#These really should be functions ...
for E in Energy:
    theta_200.append(theta_hkl(2,0,0))
    theta_222.append(theta_hkl(2,2,2))
    theta_111.append(theta_hkl(1,1,1))

#600 error with the math domain - goes over 90 deg below 5.8 keV
Energy600 = np.linspace(5.8, 10, 1000)*1000
for E in Energy600:
    theta_600.append(theta_hkl(6,0,0))

# Data tip and axis i.e., fig, ax plt.subplots(1), markers format string
plt.figure(1)
plt.plot(Energy, theta_200, label='(200)')
plt.plot(Energy, theta_222, label='(222)')
#plt.plot(Energy, theta_2_22, label='(2-22)')
plt.plot(Energy, theta_111, label='(111)')
plt.plot(Energy, theta_600, label='(600)')
#plt.title('title')
plt.xlabel('E (eV)')
plt.ylabel(r'$\theta$ (deg)')
plt.xlim(3.56*1000,10*1000)
plt.title("Bragg's law")
plt.legend()
plt.grid()
plt.savefig('figures/pre_xrd/theta(E,h,k,l).png', format='png', dpi=300)
plt.show()

#Probably start a new script from here...

# %% Scattering vectors for different diffraction angles - plot as 
# vertical lines to the structure factor difference
# \delta k = q = G Laue for diffraction

def q_hkl(h,k,l):
    return 2*pi/(d_hkl(h,k,l)*1E10) #In Ångströms

# q_200 = q_hkl(2,0,0)
# q_222 = q_hkl(2,2,2) # = q_2-22
# q_600 = q_hkl(6,0,0)

#q_111 = q_hkl(1,1,1)

# %% Atomic form factor for the structure factor
# Ref. http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php and reference therein (Journal of crystallography)

#Indium atom
def f_In(q):
    # Constants for the Gaussian parametrization of the atomic form factor
    a_In = np.array([19.1624, 18.5596, 4.2948, 2.0396])
    b_In = np.array([0.5476, 6.3776, 25.8499, 92.8029])
    c_In = 4.9391
    return np.sum(a_In*np.exp(-1*b_In*(q/(4*pi))**2)) + c_In

print(f_In(q_hkl(2,2,2)))

#Antimony atom 
def f_Sb(q):
    # Constants for the Gaussian parametrization of the atomic form factor
    a_Sb = np.array([19.6418, 19.0455, 5.0371, 2.6827])
    b_Sb = np.array([5.3034, 0.4607, 27.9074, 75.2825])
    c_Sb = 4.5909
    return np.sum(a_Sb*np.exp(-1*b_Sb*(q/(4*pi))**2)) + c_Sb

print(f_Sb(q_hkl(2,0,0)))
#print(f_Sb(2,0,0))

# %% Plotting atomic form factors over a q-range

f_In1 = []; f_Sb1 = []

interval = np.linspace(0,25,10000)
for q in interval:
    f_In1.append(f_In(q))
    f_Sb1.append(f_Sb(q))
    # print(f_Sb(q))
    # print(f_In(q))

f_In_arr = np.array(f_In1)
f_Sb_arr = np.array(f_Sb1)

plt.figure(2)
plt.plot(interval, f_In_arr, label="$f_{In}$")
plt.plot(interval, f_Sb_arr, label="$f_{Sb}$")
plt.vlines(q_hkl(2,0,0), 0, 52, colors='r', linestyles='--', label='$q_{200}$')
plt.vlines(q_hkl(2,2,2), 0, 52, colors='g', linestyles='--', label='$q_{222}$')
plt.vlines(q_hkl(6,0,0), 0, 52, colors='b', linestyles='--', label='$q_{600}$')
#r'$\theta$ (deg)'
#plt.title('title')
plt.xlabel('$q$ (Å$^{-1}$)')
#r'$\theta$ (deg)'
plt.ylabel('$f(q)$')
plt.legend()
plt.xlim(0, 25)
plt.ylim(0, 52)
plt.title("Atomic form factors over q-range")
#plt.ylim(0, max(f_In_arr))
plt.grid()
plt.savefig('figures/pre_xrd/form_factors.png', format='png', dpi=300)
plt.show()


"""Probably remove"""
#plt.figure(3)
# fig, ax = plt.subplots()
# #ax.tick_params(direction='in')
# ax.plot(interval, f_Sb_arr - f_In_arr)
# #Secondary axis
# secax = ax.secondary_xaxis('top')
# secax.set_xticks([q_hkl(2,0,0), q_hkl(2,2,2), q_hkl(6,0,0)], ['1','2','3'])
# ax.grid()
# plt.show()
"""Yep"""

"""Difference between the form factors
# The structure factor is of the form: F_hkl = 4(f_In - f_Sb) for reflections: h + k + l = 4N + 2 where N is an integer"""

plt.figure(4)
plt.plot(interval, f_Sb_arr - f_In_arr)
plt.vlines(q_hkl(2,0,0), 0, 2, colors='r', linestyles='--', label='$q_{200}$')
plt.vlines(q_hkl(2,2,2), 0, 2, colors='g', linestyles='--', label='$q_{222}$')
plt.vlines(q_hkl(6,0,0), 0, 2, colors='b', linestyles='--', label='$q_{600}$')
plt.legend()
plt.xlim(0, 25)
plt.ylim(0, f_Sb_arr[0] - f_In_arr[0])
plt.grid()
plt.savefig('figures/pre_xrd/delta_form_factors.png', format='png', dpi=300)
plt.show()

# %% Construct the crystal structure and visualize
a = 6.479E-10 # Lattice constant [m]

# Note the atomic position should be multiplied by lattice const. a
# 4 Indium atoms in an fcc crystal
x_In = a*np.array([0, 0, 1/2, 1/2])
y_In = a*np.array([0, 1/2, 0, 1/2])
z_In = a*np.array([0, 1/2, 1/2, 0])

# 4 Antimony atoms in an fcc crystal
x_Sb = x_In + a*1/4
y_Sb = y_In + a*1/4
z_Sb = z_In + a*1/4

# 3D plot of the crystal structure

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(0,0,0)
ax.scatter(x_In, y_In, z_In, s=50)
ax.scatter(x_Sb, y_Sb, z_Sb, s=50) #, c='k')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#The default axis is a bit weird but keep it like that anyways
# Can flip for example ax.set_ylim3d(5e-10,0)
ax.set_xlim3d(0, 5e-10)
ax.set_ylim3d(0, 5e-10)
ax.set_zlim3d(0, 5e-10)
#ax.grid(False)
plt.savefig('figures/pre_xrd/primitive_crystal.png', format='png', dpi=300)
plt.show()


# %% Construct the crystal structure from the translation vector - Move to a new file and integrate with the xrd-diff
"""Skip this for now, but create a function of this that takes n1, n2 and n3 as inputs and import it to the trxd.py file"""
# Note: Multiply everything by the lattice constant
pos_In1 = np.array([0,0,0])
pos_Sb1 = pos_In1 + a*1/4

#primitive translations
a1 = a*1/2*np.array([0,1,1])
a2 = a*1/2*np.array([1,0,1])
a3 = a*1/2*np.array([1,1,0])

# Evaluates the inner loop first
pos_In_list = []; pos_Sb_list = []

#Creates 2(n3 * n2 * n1) atoms
# 222 plots fcc primitive cell of In and Sb
# Note: Cannot exactly create the same 4 atom position configuration as in 
# the above cell
for n3 in range(2):
    for n2 in range(2):
        for n1 in range(2):
            pos_In_list.append(pos_In1 + (n1*a1 + n2*a2 + n3*a3))
            pos_Sb_list.append(pos_Sb1  + (n1*a1 + n2*a2 + n3*a3))
            # print('n1 = ', n1, 'n2 = ', n2, 'n3 = ', n3)
            # t = (n1*a1 + n2*a2 + n3*a3)

# This conversion is stupid, should directly work with numpy arrays
# maybe use np.empty()
# empty_array = np.array([])
# to_append = np.array([1,2,3])
# combined_array = np.append(empty_array,to_append)

array_In = np.array(pos_In_list)
array_Sb = np.array(pos_Sb_list)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x,y,z)
ax.scatter(array_In[:,0], array_In[:,1], array_In[:,2], s=50)
ax.scatter(array_Sb[:,0], array_Sb[:,1], array_Sb[:,2], s=50) 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.savefig('figures/pre_xrd/arbitrary_crystal.png', format='png', dpi=300)
#ax.grid(False)
plt.show()

# %%
