# %%
# Plot the lifetime from Mahan eq. 7.47 as a function of temperature
# This described the absorption of a LO phonon by an electron at 0 momentum (0 electron energy)
# I need to also consider that this electron re-emits the same phonon after some time
#   - This comes from considering the emission term in Fermi's golden rule
#       - Can use Andreas Wackers solution but need to consider also the Bose-Einstein distribution
# In the end, the electrons and phonons reach a thermal equilibrium and the number of phonons is then constant i.e., same as the number of phonons prior to the laser pulse?.
#   - How would I ever extract the damping constant from here? or are we rather interested how long after the laser pulse the thermal equilibrium is reached (for which the 222 bragg reflection intensity is low)?
# This method would also give me the damping constant if I fit the initial 
# xrd intensity and the final xrd intensity at some time. 

# Also a electron of energy E = hbar*omega_0 can absorb another phonon
# Note Mahan eq. 7.36 gives the number of phonons
import numpy as np
import matplotlib.pyplot as plt

hbar = 6.582119569E-16 
alpha = 0.023
k_b = 8.617333262145E-5 # eV/K

alpha = 0.0222 # Fröhlich coupling constant
omega_0 = 2*np.pi*5.76E12 # Angular frequency of a LO phonon

#%%
def lifetime(T): #This is the phonon lifetime
    a = hbar*omega_0/(k_b*T)
    #a = 1
    N_0 = 1/(np.exp(a)-1)
    return 1/(2*alpha*omega_0*N_0)


import matplotlib.pyplot as plt 
T = np.linspace(90,500,5000)

plt.figure(1)
plt.xlim(90, max(T))
plt.ylim(0,max(lifetime(T)))
plt.plot(T,lifetime(T))
plt.xlabel("T (K)")
plt.ylabel("Lifetime of LO phonon (s)")
plt.show()
#%% Use the Eq. 7.46 in Mahan to see how the absorption rate varies if the 
# electron has already absorbed a phonon

def lifetime1(T,e_p):
    a = hbar*omega_0/(k_b*T)
    N_0 = 1/(np.exp(a)-1)
    h = np.log((np.sqrt(omega_0 + e_p) + np.sqrt(e_p))/((np.sqrt(omega_0 + e_p) - np.sqrt(e_p))))
    Gamma = alpha*N_0*np.sqrt(omega_0**3/e_p)*h
    return 1/Gamma

plt.figure(2)
plt.xlim(90, max(T))
plt.ylim(0, max(lifetime1(T,7.5*omega_0)+1E-12))

#plt.plot(T,lifetime(T), label="Mahan Eq. 7.47")
plt.plot(T,lifetime1(T,0.000001*omega_0),label="e_p = 0")
plt.plot(T,lifetime1(T,1*omega_0), label="e_p = hbar*omega_0")
plt.plot(T,lifetime1(T,2*omega_0), label="e_p = 2*hbar*omega_0")
plt.plot(T,lifetime1(T,7.5*omega_0), label="e_p = 7.5*hbar*omega_0 \n = 0.18 eV (E_gap)")

plt.xlabel("T (K)")
plt.ylabel("Lifetime of LO phonon (s)")
plt.legend()
plt.show()
# %% Derive an expression for the emission rate wacker 2.11
# I think I need to consider the absorption term in fermis-golden rule for electron of 0 energy and then also consider the emission term for an electron with energy of hbar*omega_0 i.e., absorbed one phonon (electron of 0 energy can't emit a phonon) and somehow these two rates give the thermal equilibrium 

# But where do these electrons scatter as there are no states at Gamma(k=0) with these energies...

#This should be the same as my solution but with factor N_0 + 1 infront.

#%% Fundamental Constants 
echarge = 1.60217662E-19 # Coulomb
#hbar = 6.582119569E-16 #eVs
m_e = 9.10938356E-31 # kg
c = 3E10 # cm/s
eps_0 = 8.8541878128E-12 #F/m

#%% InSb material parameters
# Taken from Marder, but should replace with my own DFPT simulations
    # - ph_2.out is my most accurate simulation atm
    # - need to check LST relation omega_LO^2/omega_TO^2 = eps(0)/eps(inf)

E_gap = 0.18 # eV at 300 K
m_eff = 0.014 * m_e

eps_inf = 15.68
eps_static = 17.88
#omega = 2*pi*f
# w_TO = 2*np.pi*c*185 # rad/s
# w_LO = 2*np.pi*c*197

# GaAs material parameters - To reference
# eps_inf = 10.92
# eps_static = 12.9
# m_eff = 0.067*m_e
# E_gap = 1.44
# #E_gap = 0.06
# w_LO = 5.469E13 # Rad/s
#These values should give me 4.5E12 -> 220 fs e-ph time

#Need to use the same alphas and frequencies
alpha_check = (echarge**2/hbar)*np.sqrt(0.014*m_e/(2*hbar*omega_0))*(1/eps_inf - 1/eps_static)
alpha_marder = 1.44E8*(1/eps_inf - 1/eps_static)*np.sqrt(0.014/omega_0)
print(alpha_check, alpha_marder)

#%% Fröhlich coupling between electrons and longitudinal optical phonons
    # - Note this assumes low temperature s.t. no phonons are thermally excited
    # - I.e., neglecting the stimulated emission by electrons in Fermi's golden rule
    # - This is not completely valid as energy of a LO phonon in InSb is 24 meV and k_BT(300K) = 25 meV
    # - See Mahan book for better theories

#Note this coupling rate is twice what I have in my H.O. eqn - remove the factor of two
beta = hbar * omega_0
alpha = hbar**2/(2*m_eff)
k = np.sqrt(2*m_eff*E_gap/(hbar**2))
k_0 = np.sqrt(k**2 - (beta)/alpha)

a = echarge**2/(hbar*(2*np.pi)**2)
b = beta/(2*eps_0)*(1/eps_inf - 1/eps_static)
d = np.pi/(alpha*k)
#f = np.log(abs((k+k_0)/(k-k_0)))
f = np.log(k+k_0)-np.log(k-k_0)
Gamma = a*b*d*f

#Gamma = echarge**2/(hbar*(2*np.pi)**2)*((beta)/(2*eps_0)*(1/eps_inf - 1/eps_static))*np.pi/(alpha*k)*np.log(abs((k+k_0)/(k-k_0)))

print(f"Gamma = {Gamma} 1/s, 1/Gamma = {1/Gamma} s")