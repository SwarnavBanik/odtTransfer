"""
###############################################################################
Code for evaluating Optical Dipole Trap Potential energy
###############################################################################
Created:    Swarnav Banik  on  May 12, 2023 
"""
import numpy as np
# import src.frmtFig as frmtFig 
# import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#clrPts = frmtFig.frmtFig()

# %% Fundamental Constants ####################################################
hbar = 1.05457182 * 1e-34               # hbar [m^2 kg/ s]
kB   = 1.380649 * 1e-23                 # Boltzmann Constant [m^2 kg s^-2 K^-1]
c    = 299792458                        # Speed of sound in vacuum [m/s]
g    = 9.8                              # Acceleration due to gravity [m/s^2]

# %% Common ODT functions #####################################################
# ODT depth
def odtU_depth(m, nu = None, r = None, Gamma = None, I = None, omega0 = None, omegaL = None):
    if Gamma != None and I != None and omega0 != None and omegaL != None:
        U0 = 3*np.pi*c**2/(2*omega0**3) * (Gamma/(omega0-omegaL) + Gamma/(omega0+omegaL)) * I
    elif nu != None and r != None:
        if np.shape(nu) != np.shape(r):
            raise Exception('odt:odtU_depth: nu and r need to be of the same length.')
        U0 = 0
        for jj in range(len(r)):
            U0 = U0 + 0.5 * m * (2*np.pi*nu[jj])**2 * r[jj]**2
    else:
        raise Exception('odt:odtU_depth: please provide some inputs to process the trap depth.')
    return abs(U0)

def freqVSpower(nu0, P0, P):
    return nu0 * np.sqrt(P/P0)


# %% ODT Potentials ###########################################################
# ODT Potential: gravity
def gravityU(z,m):
    return m*g*z

# ODT Potential: Harmonic Oscillator
def odtU_HO(m, nu_x, nu_y, nu_z, x0, y0, z0_red, x,y,z):
    V = + 0.5 * m * (2*np.pi*nu_x)**2 * (x-x0)**2\
        + 0.5 * m * (2*np.pi*nu_y)**2 * (y-y0)**2\
        + 0.5 * m * (2*np.pi*nu_z)**2 * (z-z0_red)**2
    return V

# ODT Potential: Box Trap
def odtU_box(U0, x0, y0, z0_red, Lx, Ly, Lz, x,y,z):
    V = np.zeros(np.shape(x))
    indices = np.where(np.logical_or(\
                       np.logical_or(abs(x-x0) >= Lx/2, abs(y-y0) >= Ly/2),\
                                     abs(z-z0_red) >= Lz/2))
    V[indices] = U0
    return V

# ODT Potential: Disk Trap
def odtU_disk(U0, x0, y0, z0_red, Lxy, Lz, x,y,z):
    V = np.zeros(np.shape(x))
    indices = np.where(np.logical_or(np.sqrt( (x-x0)**2+(y-y0)**2 ) >= Lxy/2, abs(z-z0_red) >= Lz/2))
    V[indices] = U0
    return V

# ODT Potential: quadrupole magnetic trap
def odtU_quadMagnetic(Bdash, fieldSplit, x,y,z):
    if np.size(Bdash) != 1:
        raise Exception('odt:: odtU_quadMagnetic: bdash needs to be of size 1')
    r = 0
    if not x == []:
        r = r + (x**2)/4
    if not y == []:
        r = r + (y**2)/4
    if not z == []:
        r = r + z**2
    r = np.sqrt(r)
    #r = np.sqrt((x**2)/4 + (y**2)/4 + z**2)    
    return +fieldSplit*(Bdash*100)*r*hbar*2*np.pi*1e6

# ODT Potential: red-detuned Single beam
def odtU_redDet_singleBeam(U0, w0, wl, x0,y0,z0, x,y,z):
    zR = np.pi*w0**2/wl
    return -U0 * np.exp( -2*((y-y0)**2 + (z-z0)**2)/w0**2 ) #* np.exp( -2*((x-x0)**2 )/zR**2 )

# ODT Potential: blue-detuned ring
def odtU_blueDet_ring(U0, w0, z0_red, x,y,z):
    r = np.sqrt(x**2+y**2)
    return U0 * np.exp( -2*(r**2 + (z-z0_red)**2)/w0**2 )

# %% Test Script ##############################################################
# import na as na
# import frmtFig as frmtFig 
# import matplotlib.pyplot as plt
# from matplotlib.gridspec import GridSpec
# clrPts = frmtFig.frmtFig()

# # Inputs ======================================================================
# N_grid   = 201              # Spatial Grid Size [odd int]
# N_temp   = 200              # Temperature grid size [int]
# T_start  = 100e-9           # Start of the temperature range to be explored [K]
# T_end    = 1000e-6          # End of the temperature range to be explored [K]
# P0_red         = 50e-6      # Total P_red [W]
# trapFrqYZ0_red = 20         # Trap Frequency [Hz]
# odtWL_red      = 1064e-9    # Wavelength [m]
# w0_red         = 25*1e-6    # Beam waist size [m]
# z0_red         = 1*w0_red   # Beam waist vertical position [m]
# Bdash    = 260              # Field Gradient [G/cm]
# # Initialize ==================================================================
# x = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# y = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# z = np.linspace(-10*w0_red, +10*w0_red, N_grid)     # Spatial Grid
# X, Y, Z = np.meshgrid(x, y, z, indexing = 'ij')     # Spatial Grid
# # Compute =====================================================================
# P_red = 2
# trapFrqYZ_red = freqVSpower(trapFrqYZ0_red, P0_red, P_red)
# U0  = odtU_depth(na.m, nu = [trapFrqYZ_red, trapFrqYZ_red*odtWL_red/np.pi/w0_red], \
#         r = [w0_red, np.pi*w0_red**2/odtWL_red])
# print(f'Trap Depth of Single IR beam (f = {trapFrqYZ_red:.0f} Hz, w0 = {w0_red*1e6:.1f} um) = {U0*1e6/kB:.4f} uK')
# U_red = odtU_redDet_singleBeam(U0, w0_red, 0,0,z0_red, X, Y, Z)
# U   = + U_red \
#       + odtU_quadMagnetic(Bdash, na.fieldSplit, X,Y,Z) \
#       + gravityU(Z, na.m)
# # Plot ========================================================================  
# plot_min = np.min(U)*1e6/kB - 0.1*( np.max(U)-np.min(U) )*1e6/kB
# plot_max = np.max(U)*1e6/kB + 0.1*( np.max(U)-np.min(U) )*1e6/kB
    
# fig = plt.figure(1, figsize = (12*2,4*6))
# gs = GridSpec(2,2)
# axs = fig.add_subplot(gs[0,0:2]) 
# idx_x = np.abs(X[:,0,0]- 0).argmin()
# idx_y = np.abs(Y[0,:,0]- 0).argmin()
# axs.plot(Z[idx_x,idx_y,:]*1e6, U_red[idx_x,idx_y,:]*1e6/kB,'--', linewidth = 3, \
#           label = 'U$_{\\rm rd}$'+'(x = {:.0f} $\\mu$m, y = {:.0f} $\\mu$m)'.format(X[idx_x,0,0]*1e6,Y[0,idx_y,0]*1e6),\
#           color = clrPts[1], markersize = 10)
# axs.plot(Z[idx_x,idx_y,:]*1e6, U[idx_x,idx_y,:]*1e6/kB,'-', linewidth = 3,\
#           label = 'U$_{\\rm total}$'+'(x = {:.0f} $\\mu$m, y = {:.0f} $\\mu$m)'.format(X[idx_x,0,0]*1e6,Y[0,idx_y,0]*1e6),\
#           color = clrPts[0], markersize = 10)
# axs.axvline(z0_red*1e6, ls ='--', color = (1,1,1), linewidth = 4, label = 'IR beam center')
# axs.grid('major')
# axs.set_xlabel('Z ($\\mu$m)', fontsize = 22)
# axs.set_ylabel('(U - E$_0$)/k$_{\\rm B}$ ($\\mu$K)', fontsize = 22)
# h1, l1 = axs.get_legend_handles_labels()
# leg = axs.legend(h1, l1, loc = 'lower right')
# axs.set(ylim = (plot_min,plot_max))
# #axs.set(ylim = (0,1000))

# axs = fig.add_subplot(gs[1,0]) 
# idx = np.abs(X[:,0,0]- 0).argmin()
# axs.set_aspect('equal')
# c = axs.pcolor(Y[idx,:,0]*1e6, Z[idx,0,:]*1e6, np.transpose(U[idx,:,:])*1e6/kB, cmap='Oranges', \
#                 vmin = plot_min, vmax=plot_max)
# axs.grid('major')
# fig.colorbar(c, ax = axs, label = '$U/k_{\\rm B}$ ($\\mu$K)',orientation='horizontal')
# axs.set_ylabel('z ($\\mu$m)', fontsize = 22)
# axs.set_xlabel('y ($\\mu$m)', fontsize = 22)
# axs.set_title('U(y,z)/k$_{\\rm B}$ ($\\mu$K) at x = 0', fontsize = 22)

# axs = fig.add_subplot(gs[1,1]) 
# idx = np.abs(Z[0,0,:]- 0).argmin()
# axs.set_aspect('equal')
# c = axs.pcolor(X[:,0,idx]*1e6, Y[0,:,idx]*1e6, U[:,:,idx]*1e6/kB, cmap='Oranges', \
#                 vmin = plot_min, vmax=plot_max)
# axs.grid('major')
# fig.colorbar(c, ax = axs, label = '$U/k_{\\rm B}$ ($\\mu$K)',orientation='horizontal')
# axs.set_ylabel('x ($\\mu$m)', fontsize = 22)
# axs.set_xlabel('y ($\\mu$m)', fontsize = 22)
# axs.set_title('U(x,y)/k$_{\\rm B}$ ($\\mu$K) at z = 0', fontsize = 22)

