#%%
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from pre_xrd import q_hkl, f_In, f_Sb

#%% Build the phonon mode as a result of the (2,-2,2) laser field
# Some problem with the scaling again in the lattice constat...
a = 6.479E-10 # Lattice constant [m]

#In units of bond length
#a = 1
# scale the displacement
#a_scaling = 6.479E-10

"Equilibrium positions"
# Note the atomic position should be multiplied by lattice const. a
# 4 Indium atoms in an fcc crystal
x_In = np.array([0, 0, 1/2, 1/2])
y_In = np.array([0, 1/2, 0, 1/2])
z_In = np.array([0, 1/2, 1/2, 0])

# 4 Antimony atoms in an fcc crystal
x_Sb = x_In + 1/4
y_Sb = y_In + 1/4
z_Sb = z_In + 1/4

"Simulated atomic motion from harmonic oscillator equations as a response to an electric field in ho_solve_odeint.py file"
"Total field strength 100 kV/cm with polarization along (2,-2,2) direction"

# For Indium atoms in the crystal
# Datasets contain 10000 points over a time interval of 50 ps
dirx_In = np.loadtxt("data/dirx.txt")/a #Divide by a to get in units of the lattice const
diry_In = np.loadtxt("data/diry.txt")/a #Does the lattice constant cancel in the exponent
dirz_In = np.loadtxt("data/dirz.txt")/a #Of the structure factor???

#If i divide by the lattice constant, is the motion still correct?

# The movement of Antimony atoms is in the opposite direction scaled by the mass ratio
# (Born effective charge and Raman tensors are opposite signs for In and Sb)
# Also this is the optical T-mode
# Scaling factor
M_In = 114.82; M_Sb = 121.71 # amu
C = M_In/M_Sb

dirx_Sb = C * (-dirx_In)
diry_Sb = C * (-diry_In)
dirz_Sb = C * (-dirz_In)

# Create vectors containing the atomic positions in the crystal i.e., 
# mode_x_In = x_In + dirx_In 
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

    #print(mode_x_In)

# Convert to numpy arrays
mode_x_In_arr = np.array(mode_x_In) #mode_x_In[0,:] is In atom one 
mode_y_In_arr = np.array(mode_y_In)
mode_z_In_arr = np.array(mode_z_In)

mode_x_Sb_arr = np.array(mode_x_Sb)
mode_y_Sb_arr = np.array(mode_y_Sb)
mode_z_Sb_arr = np.array(mode_z_Sb)


#%% TIME-RESOLVED XRD INTENSITY
"""create a structure factor function that takes h,k,l as an input"""
# The sum is over the 4 atoms, but the structure factors should return an array
# of length 10000

#Need to fix the ranges to agree with steps in ho_solve_odeint.py
def struct_F_In123(h,k,l):
    f_In123 = []
    for i in range(10000): 
        f_In123.append(f_In(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_In_arr[:, i])
                                                               + k*(mode_y_In_arr[:, i])
                                                               + l*(mode_z_In_arr[:, i])))))

    return np.array(f_In123)

#maybe need to convert the lists to an array

def struct_F_Sb123(h,k,l):
    f_Sb123 = []
    for i in range(10000): 
        f_Sb123.append(f_Sb(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_Sb_arr[:, i])
                                                               + k*(mode_y_Sb_arr[:, i])
                                                               + l*(mode_z_Sb_arr[:, i])))))

    return np.array(f_Sb123)


def struct_F(h,k,l):
    return struct_F_In123(h,k,l) + struct_F_Sb123(h,k,l)

def Intensity(h,k,l):
    return struct_F(h,k,l) * np.conj(struct_F(h,k,l))

t_int = np.linspace(0,10,10000)*1E-12

plt.figure()
plt.plot(t_int, Intensity(2,2,2), label="(2,2,2)")
plt.plot(t_int, Intensity(2,-2,2), label="(2,-2,2)")
plt.plot(t_int, Intensity(2,0,0), label="(2,0,0)")
plt.plot(t_int, Intensity(6,0,0), label="(6,0,0)")
plt.xlabel("t (s)")
plt.ylabel("Intensity (arb. units)")
plt.xlim(0, max(t_int))
plt.legend()
plt.grid(True)
#plt.savefig('figures/trxd/TRXD.png', format='png', dpi=300)
plt.show()

# print(struct_F(2,2,2))
# print(Intensity(2,2,2))

# %% Plots for proposal

plt.figure()
#Need to start plotting 
plt.plot(np.linspace(0.6E-12,10E-12,9400), Intensity(2,-2,2)[600:])
#plt.plot(t_int, Intensity(2,-2,2))
plt.xlabel("t (ps)")
plt.ylabel("Intensity (arb. u.)")
plt.rcParams.update({'font.size': 20})

plt.xlim(0.6E-12, 10E-12)
plt.xticks([0.2E-11, 0.4E-11, 0.6E-11, 0.8E-11, max(t_int)], (2, 4, 6, 8, 10)) #[pm]
#For 20 ps range
#plt.xticks([0.5E-11, 1.0E-11, 1.5E-11, max(ts)])
#plt.xticks(np.arange(0.5E-12,2E-11,0.5E-11)) 

# plt.ylim(-1.6E-13, 1.6E-13)
# #plt.yticks([-1.0E-13, -0E-13, 1.0E-13], (-0.1, 0, 0.1)) #[pm]
# plt.yticks([-1.5E-13, -1.0E-13, -0.5E-13, -0E-13, 0.5E-13, 1.0E-13, 1.5E-13], (-0.15, '0.10', 0.05, 0, 0.05, '0.10', 0.15))

#For 20 ps plot
#plt.yticks(np.arange(-1.5E-13, 1.55E-13, 0.5E-13))

#plt.legend()
plt.grid(True)
#plt.savefig('figures/trxd/TRXD.png', format='png', dpi=300)
plt.show()

# %%
