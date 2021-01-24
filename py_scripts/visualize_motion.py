#%% Loading packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# Multiply the movement by 50 for visualization purposes
dirx_In = 50*np.loadtxt("data/dirx.txt") 
diry_In = 50*np.loadtxt("data/diry.txt") 
dirz_In = 50*np.loadtxt("data/dirz.txt") 

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
# DIVIDE BY a
mode_x_In_arr = np.array(mode_x_In) #mode_x_In[0,:] is atom one 
mode_y_In_arr = np.array(mode_y_In)
mode_z_In_arr = np.array(mode_z_In)

mode_x_Sb_arr = np.array(mode_x_Sb)
mode_y_Sb_arr = np.array(mode_y_Sb)
mode_z_Sb_arr = np.array(mode_z_Sb)

#%% # (and the decay in 50 ps time interval as a response to the laser field)

# Don't want to plot all 10000 images for the gif.
xx = 20
no_images = int(10000/xx)
# Multiplier for the atomic motion (otherwise dont see anything)
# mult = 10000 # This just multiplies the values but not the motion ... 
# Need to multiply the dirx_In, diry_In, dirz_In

for i in range(no_images):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(mode_x_In_arr[:, xx*i], mode_y_In_arr[:, xx*i], mode_z_In_arr[:, xx*i], s=50)
    ax.scatter(mode_x_Sb_arr[:, xx*i], mode_y_Sb_arr[:, xx*i], mode_z_Sb_arr[:, xx*i], s=50)
    # The axis is a bit weird ... 
    # ax.set_xlim3d(0, 1)
    # ax.set_ylim3d(0, 1) 
    # ax.set_zlim3d(0, 1)
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
    #print(filename)
    images.append(imageio.imread(filename))
# Also the axis are quite weird
# Also could import the translation vector crystal structure instead, i.e., more atoms in the cell
imageio.mimsave('figures/trxd/gif/Motion_multiplied_by_50.gif', images)

#Can I plot a time from 0 to 50 ps in 10000/xx steps to the image
#%% Static visualization - delete maybe (probably useless)

# # Probably should delete as I have the gif in next cell.
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# #ax.scatter(0,0,0)
# ax.scatter(x_In, y_In, z_In, s=50)
# # Problem: How to add the 10000 datapoints to the 4 entries of x_In? need to do the same thing in the structure factors also
# # Probably create a new vector already above
# # Matlab understands the below, 
# ax.scatter(x_In, y_In, z_In, s=50)
# #This can be used to visualize the movement
# ax.scatter(mode_x_In_arr[:, 100], mode_y_In_arr[:, 100], mode_z_In_arr[:, 100], s=30, c='k')

# ax.scatter(x_Sb, y_Sb, z_Sb, s=50)
# #This can be used to visualize the movement
# ax.scatter(mode_x_Sb_arr[:, 100], mode_y_Sb_arr[:, 100], mode_z_Sb_arr[:, 100], s=30, c='k')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# #ax.grid(False)
# plt.savefig('figures/trxd/phonon_mode_crystal.png', format='png', dpi=300)
# plt.show()