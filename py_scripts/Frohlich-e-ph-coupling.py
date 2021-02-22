# Need to make this work, but the problem is that this is defined for a k-point and 
# in my case we are probing gamma i.e. k=0
#%%
import numpy as np
#from math import pi, log, sqrt, exp

#%% Fundamental Constants 
echarge = 1.60217662E-19 # Coulomb
hbar = 6.582119569E-16 #eVs
m_e = 9.10938356E-31 # kg
c = 3E10 # cm/s
eps_0 = 8.8541878128E-12 #F/m

#%% InSb material parameters
# Taken from Marder, but should replace with my own DFPT simulations
    # - ph_2.out is my most accurate simulation atm
    # - need to check LST relation omega_LO^2/omega_TO^2 = eps(0)/eps(inf)

E_gap = 0.18 # eV at 300 K
m_eff = 0.014*m_e

eps_inf = 15.68
eps_static = 17.88
#omega = 2*pi*f
w_TO = 2*np.pi*c*185 # rad/s
w_LO = 2*np.pi*c*197

#%% GaAs material parameters - To reference
# eps_inf = 10.92
# eps_static = 12.9
# m_eff = 0.067*m_e
# E_gap = 1.44
# #E_gap = 0.06
# w_LO = 5.469E13 # Rad/s
#These values should give me 4.5E12 -> 220 fs e-ph time

#%% Fröhlich coupling between electrons and longitudinal optical phonons
    # - Note this assumes low temperature s.t. no phonons are thermally excited
    # - I.e., neglecting the stimulated emission by electrons in Fermi's golden rule
    # - This is not completely valid as energy of a LO phonon in InSb is 24 meV and k_BT(300K) = 25 meV
    # - See Mahan book for better theories

#Note this coupling rate is twice what I have in my H.O. eqn - remove the factor of two
beta = hbar * w_LO
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
# %%
