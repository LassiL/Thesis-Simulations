#%%
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from pre_xrd import q_hkl, f_In, f_Sb

#%% Build the phonon mode as a result of the (2,-2,2) laser field
a = 6.479E-10 # Lattice constant [m]

"Equilibrium positions"
# Note the atomic position should be multiplied by lattice const. a
# 4 Indium atoms in an fcc crystal
x_In = a*np.array([0, 0, 1/2, 1/2])
y_In = a*np.array([0, 1/2, 0, 1/2])
z_In = a*np.array([0, 1/2, 1/2, 0])

# 4 Antimony atoms in an fcc crystal
x_Sb = x_In + a*1/4
y_Sb = y_In + a*1/4
z_Sb = z_In + a*1/4

"Simulated atomic motion from harmonic oscillator equations as a response to an electric field in ho_solve_odeint.py file"
"Total field strength 100 kV/cm with polarization along (2,-2,2) direction"

# For Indium atoms in the crystal
# Datasets contain 10000 points over a time interval of 50 ps
dirx_In = 100*np.loadtxt("data/dirx.txt")
diry_In = 100*np.loadtxt("data/diry.txt")
dirz_In = 100*np.loadtxt("data/dirz.txt")

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
mode_x_In_arr = np.array(mode_x_In) #mode_x_In[0,:] is atom one 
mode_y_In_arr = np.array(mode_y_In)
mode_z_In_arr = np.array(mode_z_In)

mode_x_Sb_arr = np.array(mode_x_Sb)
mode_y_Sb_arr = np.array(mode_y_Sb)
mode_z_Sb_arr = np.array(mode_z_Sb)

#%% Visualize the crystal movement in the crystal
# Not very good ... but seems like they are moving correct
# Need to make a gif and more atoms
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(0,0,0)
ax.scatter(x_In, y_In, z_In, s=50)
# Problem: How to add the 10000 datapoints to the 4 entries of x_In? need to do the same thing in the structure factors also
# Probably create a new vector already above
# Matlab understands the below, 
ax.scatter(x_In, y_In, z_In, s=50)
#This can be used to visualize the movement
ax.scatter(mode_x_In_arr[:, 100], mode_y_In_arr[:, 100], mode_z_In_arr[:, 100], s=30, c='k')

ax.scatter(x_Sb, y_Sb, z_Sb, s=50)
#This can be used to visualize the movement
ax.scatter(mode_x_Sb_arr[:, 100], mode_y_Sb_arr[:, 100], mode_z_Sb_arr[:, 100], s=30, c='k')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#ax.grid(False)
plt.savefig('figures/trxd/phonon_mode_crystal.png', format='png', dpi=300)
plt.show()

#print(f_In(q_hkl(2,0,0)))

#%% Make a gif of the motion

# Don't want to plot all 10000 images for the gif.
xx = 20
no_images = int(10000/xx)

for i in range(no_images):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(mode_x_In_arr[:, xx*i], mode_y_In_arr[:, xx*i], mode_z_In_arr[:, xx*i], s=50)
    ax.scatter(mode_x_Sb_arr[:, xx*i], mode_y_Sb_arr[:, xx*i], mode_z_Sb_arr[:, xx*i], s=50)
    # The axis is a bit weird ... 
    ax.set_xlim3d(0, 5e-10)
    ax.set_ylim3d(0, 5e-10)
    ax.set_zlim3d(0, 5e-10)
    plt.savefig('figures/trxd/gif/{}.png'.format(i), format='png', dpi=100)
    plt.close()

#Kinda sketchy to create a list like this put whatever
list_filenames = []
for i in range(no_images):
    list_filenames.append('figures/trxd/gif/{}.png'.format(i))
#print(list_filenames)

import imageio
images = []
for filename in list_filenames:
    print(filename)
    images.append(imageio.imread(filename))
#Looks good but the oscillation seems quite high??
imageio.mimsave('figures/trxd/gif/movie.gif', images)


# %% TIME-RESOLVED XRD INTENSITY
"""create a structure factor function that takes h,k,l as an input"""
# Loop the function or make a loop of the function
# Perhaps need to loop over i = 10000 in the function
# Also need to append them to an array

def struct_F_In(h,k,l):
    return f_In(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_In_arr[:, i])
                                                   + k*(mode_y_In_arr[:, i])
                                                   + l*(mode_z_In_arr[:, i]))))

def struct_F_Sb(h,k,l):
    return f_Sb(q_hkl(h,k,l))*np.sum(np.exp(2*pi*1j*(h*(mode_x_Sb_arr[:, i])
                                                   + k*(mode_y_Sb_arr[:, i])
                                                   + l*(mode_z_Sb_arr[:, i]))))

def Intensity(h,k,l):
    return (struct_F_In(h,k,l) + struct_F_In(h,k,l)) * np.conj(struct_F_In(h,k,l) + struct_F_Sb(h,k,l))
