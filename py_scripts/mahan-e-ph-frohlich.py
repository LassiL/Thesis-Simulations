
#%% Mahan eq. 7.46
import numpy as np

hbar = 6.582119569E-16 # eV*s

#%% InSb
#alpha = 0.023 #taken from marder -> should do with my own simulations
#See marder 22.41 or mahan 7.2 for calculation
omega_0 = 2*np.pi*5.76E12 #rad/s
#e_p doesnt change anything???
e_p = 0.18 # Energy of a particle at VBM (gamma point of InSb) in eV

#%% GaAs: gives 1.86E12 1/s i.e., 500 fs lifetime 
# alpha = 0.068
# omega_0 = 2*np.pi*8.5E12 #rad/s
# #e_p doesnt change anything???
# e_p = 1.44 # Energy of a particle at VBM (gamma point of GaAs) in eV

#%%
k_b = 8.617333262145E-5 #eV/K
T = 293 # K the lifetime varies very rapidly with temperature
beta = 1/(k_b*T) #approx 1/0.025 1/eV at 300 K
N_0 = 1/(np.exp(beta*hbar*omega_0)-1) #Have to define this with hbar

#How does this not depend on the energy? also im adding omega_0 (rad/s) with eV...
Gamma = alpha*N_0*np.sqrt(omega_0**3/e_p)*np.log((np.sqrt(omega_0 + e_p) + np.sqrt(e_p))/((np.sqrt(omega_0 + e_p) - np.sqrt(e_p))))

# Note that these gammas are twice my trxd gamma
print(f"Gamma = {Gamma} 1/s, 1/Gamma = {1/Gamma} s.")


# %% At zero momentum the life time is evaluated as mahan eq. 7.47
#have to divide by hbar here
Gamma_0 = 2*alpha*omega_0*N_0

print(f"Gamma_0 = {Gamma_0} 1/s, 1/Gamma_0 = {1/Gamma_0} s.")

# gives that the lifetime at 0 momentum is 1 ps at 300 K
# at 100 K the lifetime is 10 ps
# Lockwood et al. Solid State Communications 136 (2005) 404–409
# approximates the LO lifetime as 10 ps from the dielectric function ... 
# # %%
