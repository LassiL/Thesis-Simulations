#%%
"""This scripts takes the solutions of the harmonic oscillator equations given by the script ho_solve_odeint.py 
and calculates and plots the TXRD intenity of the different weak reflections using the form factors from pre_xrd"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from math import pi

#Import scattering vector and form factors from pre_xrd.py
from pre_xrd import q_hkl, f_In, f_Sb 

#%% Build the crystal structure
a = 6.479E-10 # Lattice constant [m]

# Positions of the In and Sb atoms in the zincblende conventional unit cell 
# in units of the lattice constant

# 4 Indium atoms in an fcc crystal
x_In = np.array([0, 0, 1/2, 1/2])
y_In = np.array([0, 1/2, 0, 1/2])
z_In = np.array([0, 1/2, 1/2, 0])

# 4 Antimony atoms in an fcc crystal shifted
x_Sb = x_In + 1/4
y_Sb = y_In + 1/4
z_Sb = z_In + 1/4

#%% Load the harmonic oscillator solutions for the displacements from the ho_solve_odeint.py and create the THz excited phonon modes 

# For Indium atoms in the crystal
# Divide by the lattice constant due to the structure factor taking the Cartesian positions
dirx_In = np.loadtxt("data/dirx.txt")/a 
diry_In = np.loadtxt("data/diry.txt")/a 
dirz_In = np.loadtxt("data/dirz.txt")/a 

# The movement of Antimony atoms is in the opposite direction scaled by the mass ratio
# (Born effective charge and Raman tensors are equal but opposite signs for In and Sb)
# Scaling factor
M_In = 114.82; M_Sb = 121.71 # amu
C = M_In/M_Sb

# Displacement of the Sb atoms as given by the In atoms
dirx_Sb = C * (-dirx_In)
diry_Sb = C * (-diry_In)
dirz_Sb = C * (-dirz_In)

# Create lists containing the atomic positions in the crystal as a function of time i.e., 
# equilbrium positions and the displacement: mode_x_In = x_In + dirx_In  
mode_x_In = []; mode_y_In = []; mode_z_In = []
mode_x_Sb = []; mode_y_Sb = []; mode_z_Sb = []

no_atoms = 4
for i in range(no_atoms):
    mode_x_In.append(x_In[i] + dirx_In[:])
    mode_y_In.append(y_In[i] + diry_In[:])
    mode_z_In.append(z_In[i] + dirz_In[:])

    mode_x_Sb.append(x_Sb[i] + dirx_Sb[:])
    mode_y_Sb.append(y_Sb[i] + diry_Sb[:])
    mode_z_Sb.append(z_Sb[i] + dirz_Sb[:])

# Convert to numpy arrays for easier handling
mode_x_In_arr = np.array(mode_x_In) #mode_x_In[0,:] is In atom one, mode_x_In[1,:] is In atom two 
mode_y_In_arr = np.array(mode_y_In) #etc
mode_z_In_arr = np.array(mode_z_In)

mode_x_Sb_arr = np.array(mode_x_Sb)
mode_y_Sb_arr = np.array(mode_y_Sb)
mode_z_Sb_arr = np.array(mode_z_Sb)


#%% TIME-RESOLVED XRD INTENSITY

#Note: the ranges have to agree with the timestep given in ho_solve_odeint.py, i.e.,
#ts = np.linspace(0, 11E-12, 40000) in line 86

#Structure factor arising from the 4 indium atoms 
def struct_F_In123(h,k,l):
    f_In123 = []
    for i in range(40000): 
        f_In123.append(f_In(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_In_arr[:, i])
                                                               + k*(mode_y_In_arr[:, i])
                                                               + l*(mode_z_In_arr[:, i])))))

    return np.array(f_In123)

#Structure factor arising from the 4 antimony atoms 
def struct_F_Sb123(h,k,l):
    f_Sb123 = []
    for i in range(40000): 
        f_Sb123.append(f_Sb(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_Sb_arr[:, i])
                                                               + k*(mode_y_Sb_arr[:, i])
                                                               + l*(mode_z_Sb_arr[:, i])))))

    return np.array(f_Sb123)

#The total structure factor is the sum of the structure factors of the In and Sb lattice 
def struct_F(h,k,l):
    return struct_F_In123(h,k,l) + struct_F_Sb123(h,k,l)

def Intensity(h,k,l):
    return struct_F(h,k,l) * np.conj(struct_F(h,k,l))

#%% PLOT TXRD
#Note: Defining time interval for the plots, has to match the time interval of the 
#ho_solve_odeint.py i.e., ts = np.linspace(0, 11E-12, 40000) in line 86
t_int = np.linspace(0,11E-12,40000)

plt.figure()
plt.plot(t_int, Intensity(2,-2,2), label=r"$2\bar{2}2$") 
plt.plot(t_int, Intensity(2,2,2), label=r"$222$")
plt.plot(t_int, Intensity(6,0,0), label=r"$600$")
plt.plot(t_int, Intensity(2,0,0), label=r"$200$")
plt.xlabel("t (s)")
plt.ylabel("Intensity (arb. units)")
plt.xlim(0, 10.95E-12)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('figures/trxd/TRXD.png', format='png', dpi=300)
plt.show()
# %%