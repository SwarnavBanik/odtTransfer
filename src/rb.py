"""
###############################################################################
Atomic Properties
###############################################################################
Created:    Swarnav Banik  on  May 12, 2023 
"""
import numpy as np

# %% Fundamental Constants ####################################################
hbar = 1.05457182 * 1e-34               # reduced Plancks Constant [m^2 kg/ s]
kB   = 1.380649 * 1e-23                 # Boltzmann Constant [m^2 kg s^-2 K^-1]
c    = 299792458                        # Speed of sound in vacuum [m/s]

# %% Rubidium ###################################################################
m = 85.4678*1.66054*1E-27              # Mass [kg]
fieldSplit = 0.7                       # Zeeman splitting b/w adjacent magnetic sublevels [MHz/G]