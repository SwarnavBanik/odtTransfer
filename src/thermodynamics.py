"""
###############################################################################
Code for evaluating Thermodynamic Properties of atoms in an ODT
###############################################################################
Created:    Swarnav Banik  on  May 12, 2023 
"""
import numpy as np
from scipy import integrate

# %% Fundamental Constants ####################################################
hbar = 1.05457182 * 1e-34               # reduced Plancks Constant [m^2 kg/ s]
kB   = 1.380649 * 1e-23                 # Boltzmann Constant [m^2 kg s^-2 K^-1]
c    = 299792458                        # Speed of sound in vacuum [m/s]

# %% Thermodynamic Functions ##################################################

def deBroglieWL(T,m):
    lambda_dB = (2*np.pi*hbar**2/(m*kB*T))**(1/2)
    return lambda_dB

# def effVolume(U, T, r):
#     if np.ndim(U) != len(r):
#         raise Exception('thermodynamics:: Potential U and coordinate r need to be of the same dimesion.')
#     if np.size(T) == 1:
#         T = [T]
#     V = np.zeros(np.shape(T))
#     for jj in range(len(T)):
#         integrand = np.exp(-U/kB/T[jj]) 
#         for ii in range(len(r)):
#             integrand = integrate.trapz(integrand, r[ii], axis=0)
#         V[jj] = integrand
#     return V  

# def partitionFnc(V,T,m):
#     lambda_dB = deBroglieWL(T,m)
#     return V /lambda_dB**3

# def helmholtzPerAtom(V,T,m):
#     return -kB*T*np.log(partitionFnc(V,T,m))

# def entropyPerAtom(V,T,m):
#     if np.size(T) == 1:
#         raise Exception('thermodynamics:: entropyPerAtom: T needs to be an array of size > 1')
#     A = helmholtzPerAtom(V,T,m)
#     return -np.gradient(A)/np.mean(np.diff(T))

# %% Thermodynamic properties v/s Temperature #################################

class thermoODT:
    """
    Thermodynamic Properties as a function of T of an ODT
    """
    def __init__(self, U, r, m, T_start = 0.1e6, T_end = 10e-6, N_T = 100):
        if np.ndim(U) != len(r):
            raise Exception('thermoODT: Potential U and coordinate r need to be of the same dimesion.')        
        self.T = np.linspace(T_start, T_end, N_T)   # Temperature Scale [K]
        self.U = U                                  # Potential Energy [J]
        self.r = r                                  # Coordinate Space [m]
        self.m = m                                  # Atomic Mass [kg]
        self.dim = len(r)                           # Dimension of the coordinate system
        
    def effVolume(self):
        if np.size(self.T) == 1:
            self.T = [self.T]
        V = np.zeros(np.shape(self.T))
        for jj in range(len(self.T)):
            integrand = np.exp(-self.U/kB/self.T[jj]) 
            for ii in range(len(self.r)):
                integrand = integrate.trapz(integrand, self.r[ii], axis=0)
            V[jj] = integrand
        self.V = V
        
    def partitionFnc(self):
        if not hasattr(self, 'V'):
            self.effVolume()
        lambda_dB = deBroglieWL(self.T,self.m)
        self.eta = self.V /lambda_dB**self.dim
    
    def helmholtzPerAtom(self):
        if not hasattr(self, 'eta'):
            self.partitionFnc()
        self.A = -kB*self.T*np.log(self.eta)
        
    def entropyPerAtom(self):
        if np.size(self.T) == 1:
            raise Exception('thermoODT:: entropyPerAtom: T needs to be an array of size > 1')
        if not hasattr(self, 'A'):
            self.helmholtzPerAtom()
        self.S = -np.gradient(self.A)/np.mean(np.diff(self.T))
        
    def propsAtT0(self, T0):
        if T0>np.max(self.T) or T0<np.min(self.T):
            raise Exception('thermoODT:: propsAtT0: T0 needs to be in the range of T values.')
        if not hasattr(self, 'S'):
            self.entropyPerAtom()
        idx_T0  = np.abs(self.T - T0).argmin()
        T0      = self.T[idx_T0]
        VatT0   = self.V[idx_T0]
        etaatT0 = self.eta[idx_T0]
        AatT0   = self.A[idx_T0]
        SatT0   = self.S[idx_T0]
        return T0, VatT0, etaatT0, AatT0, SatT0
    
    def propsAtS0(self, S0):        
        if not hasattr(self, 'S'):
            self.entropyPerAtom()
        if S0>np.max(self.S) or S0<np.min(self.S):
            raise Exception('thermoODT:: propsAtS0: S0 needs to be in the range of S values.')
        idx_S0  = np.abs(self.S - S0).argmin()
        T0atS0  = self.T[idx_S0]
        VatS0   = self.V[idx_S0]
        etaatS0 = self.eta[idx_S0]
        AatS0   = self.A[idx_S0]
        S0      = self.S[idx_S0]
        return T0atS0, VatS0, etaatS0, AatS0, S0

    
# %% Test Script ##############################################################
# import na as na
# import odt as odt
# import frmtFig as frmtFig 
# import matplotlib.pyplot as plt
# from matplotlib.gridspec import GridSpec
# clrPts = frmtFig.frmtFig()

# # Inputs ======================================================================
# N_grid   = 201              # Spatial Grid Size [odd int]
# N_temp   = 200              # Temperature grid size [int]
# T_start  = 100e-9           # Start of the temperature range to be explored [K]
# T_end    = 1000e-6          # End of the temperature range to be explored [K]
# w0_red   = 25*1e-6          # Beam waist size [m]
# Bdash    = 200              # Field Gradient [G/cm]
# # Initialize ==================================================================
# T = np.linspace(T_start, T_end, N_temp)             # Temperature values
# #T = 1e-3
# x = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# y = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# z = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# X, Y, Z = np.meshgrid(x, y, z, indexing = 'ij')     # Spatial Grid
# # Compute =====================================================================
# U    = odt.odtU_quadMagnetic(Bdash, na.fieldSplit, X,Y,Z)
# V    = effVolume(U, T, [x,y,z])
# eta  = partitionFnc(V, T, na.m)
# A    = helmholtzPerAtom(V, T, na.m)
# S    = entropyPerAtom(V, T, na.m)
# # Plot ========================================================================
# fig = plt.figure(1, figsize = (12*1,4*6))
# gs = GridSpec(4,1)
# axs = fig.add_subplot(gs[0,0]) 
# axs.loglog(T*1e6, V*1e6,'.', linewidth = 5, color = clrPts[0], markersize = 10)
# axs.set_ylabel('V$_0$ (cm$^3$)', fontsize = 22)
# axs.grid('major')
# axs = fig.add_subplot(gs[1,0]) 
# axs.loglog(T*1e6, eta,'.', linewidth = 5, color = clrPts[0], markersize = 10)
# axs.set_ylabel('$\\xi$', fontsize = 22)
# axs.grid('major')
# axs = fig.add_subplot(gs[2,0]) 
# axs.plot(T*1e6, A,'.', linewidth = 5, color = clrPts[0], markersize = 10)
# axs.set_xlabel('T ($\\mu$K)', fontsize = 22)
# axs.set_ylabel('A/N (J)', fontsize = 22)
# axs.grid('major')
# axs = fig.add_subplot(gs[3,0]) 
# axs.loglog(T*1e6, S,'.', linewidth = 5, color = clrPts[0], markersize = 10)
# axs.set_xlabel('T ($\\mu$K)', fontsize = 22)
# axs.set_ylabel('S/N (J/K)', fontsize = 22)
# axs.grid('major')
