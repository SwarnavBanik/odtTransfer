"""
###############################################################################
Code for testing the convergence of thermodynamic properties of atoms in an ODT
###############################################################################
Created:    Swarnav Banik  on  May 26, 2023 
"""

from src.frmtFig import frmtFig 
import matplotlib as mpl
import matplotlib.pyplot as plt
clrPts, mpl, plt = frmtFig(mpl, plt, FS_title = 22, FS_tickLabel = 20, FS_axisLabel = 22)

import numpy as np
import src.na as na
import src.odt as odt
import src.gridMnp as gridMnp
import src.definePaths as paths
import src.thermodynamics as thermo
from matplotlib.gridspec import GridSpec

# %% Local Functions ##########################################################
def plot3DU(U,X,Y,Z, x0 = 0, y0 = 0, z0 = 0, figNo = 1, plot_min = [], plot_max = [],\
            plot_l = []):
    if plot_min == []:
        plot_min = np.min(U)*1e6/thermo.kB - 0.1*( np.max(U)-np.min(U) )*1e6/thermo.kB
    if plot_max == []:
        plot_max = np.max(U)*1e6/thermo.kB + 0.1*( np.max(U)-np.min(U) )*1e6/thermo.kB
    if plot_l == []:
        plot_l = np.max(abs(Z))*1e6

    fig = plt.figure(figNo, figsize = (12*1.2,12.*1.2))
    gs = GridSpec(2,2)
    axs = fig.add_subplot(gs[0,0]) 
    idx_x = np.abs(X[:,0,0]- x0).argmin()
    idx_y = np.abs(Y[0,:,0]- y0).argmin()
    axs.plot(Z[idx_x,idx_y,:]*1e6, U[idx_x,idx_y,:]*1e6/thermo.kB,'-', linewidth = 2,\
              label = '(x = {:.0f} $\\mu$m, y = {:.0f} $\\mu$m)'.format(X[idx_x,0,0]*1e6,Y[0,idx_y,0]*1e6),\
              color = clrPts[1], markersize = 10)
    axs.grid('major')
    axs.set_xlabel('z ($\\mu$m)')
    axs.set_ylabel('(U - E$_0$)/k$_{\\rm B}$ ($\\mu$K)')
    h1, l1 = axs.get_legend_handles_labels()
    axs.legend(h1, l1, loc = 'upper right')
    axs.set(ylim = (plot_min,plot_max))
    axs.set(xlim = (-plot_l, +plot_l))
    
    axs = fig.add_subplot(gs[0,1]) 
    idx_y = np.abs(Y[0,:,0]- y0).argmin()
    idx_z = np.abs(Z[0,0,:]- z0).argmin()
    axs.plot(X[:,idx_y,idx_z]*1e6, U[:,idx_y,idx_z]*1e6/thermo.kB,'-', linewidth = 2,\
              label = '(y = {:.0f} $\\mu$m, z = {:.0f} $\\mu$m)'.format(Y[0,idx_y,0]*1e6,Z[0,0,idx_z]*1e6),\
              color = clrPts[3], markersize = 10)
    axs.grid('major')
    axs.set_xlabel('x ($\\mu$m)')
    h1, l1 = axs.get_legend_handles_labels()
    axs.legend(h1, l1, loc = 'upper right')
    axs.set(ylim = (plot_min,plot_max))
    axs.set(xlim = (-plot_l, +plot_l))

    axs = fig.add_subplot(gs[1,0]) 
    idx = np.abs(X[:,0,0]- x0).argmin()
    axs.set_aspect('equal')
    c = axs.pcolor(Y[idx,:,0]*1e6, Z[idx,0,:]*1e6, np.transpose(U[idx,:,:])*1e6/thermo.kB, cmap='Oranges', \
                    vmin = plot_min, vmax=plot_max)
    axs.grid('major')
    fig.colorbar(c, ax = axs, label = 'U({:.0f}'.format(x0*1e6)+\
                 '$\\mu$m, y,z) /k$_{\\rm B}$ ($\\mu$K)', orientation='horizontal')
    axs.axvline(y0, ls = '--', linewidth = 3, color = clrPts[1])
    axs.set_ylabel('z ($\\mu$m)')
    axs.set_xlabel('y ($\\mu$m)')
    axs.set(xlim = (-plot_l, +plot_l))
    axs.set(ylim = (-plot_l, +plot_l))

    axs = fig.add_subplot(gs[1,1]) 
    idx = np.abs(Z[0,0,:]- z0).argmin()
    axs.set_aspect('equal')
    c = axs.pcolor(X[:,0,idx]*1e6, Y[0,:,idx]*1e6, U[:,:,idx]*1e6/thermo.kB, cmap='Oranges', \
                    vmin = plot_min, vmax=plot_max)
    axs.axvline(y0, ls = '--', linewidth = 3, color = clrPts[3])
    axs.grid('major')
    fig.colorbar(c, ax = axs, label = 'U(x,y,{:.0f}'.format(z0*1e6)+\
                 '$\\mu$m) /k$_{\\rm B}$ ($\\mu$K)', orientation='horizontal')
    axs.set_ylabel('x ($\\mu$m)')
    axs.set_xlabel('y ($\\mu$m)')
    axs.set(xlim = (-plot_l, +plot_l))
    axs.set(ylim = (-plot_l, +plot_l))    
   
def plotConv(figNo,X,Y,Z,Uval,T,N,z0,l,S,V):
    integrand = np.exp(-Uval/thermo.kB/T) 
    plot_min = np.min(Uval)*1e6/thermo.kB - 0.1*( np.max(Uval)-np.min(Uval) )*1e6/thermo.kB
    plot_max = np.max(Uval)*1e6/thermo.kB + 0.1*( np.max(Uval)-np.min(Uval) )*1e6/thermo.kB
    
    fig = plt.figure(figNo, figsize = (12*1,3*6))
    gs = GridSpec(3,1)
    axs = fig.add_subplot(gs[0,0]) 
    idx_x = np.abs(X[:,0,0]-0).argmin()
    idx_y = np.abs(Y[0,:,0]-0).argmin()
    axs.plot(Z[idx_x, idx_y,:]*1e6, Uval[idx_x, idx_y,:]*1e6/thermo.kB, '-', label = 'U at (0,0,z)', \
             color = clrPts[0], linewidth = 2)
    axs.set(ylim = (plot_min,plot_max))
    axs2 = axs.twinx()
    axs2.plot(Z[idx_x, idx_y,:]*1e6, integrand[idx_x, idx_y,:], '-', label = 'e$^{-U/k_{\\rm B}T}$ at (0,0,z)',\
              color = clrPts[1], linewidth = 2)
    axs2.set(ylim = (0,1))
    h1, l1 = axs.get_legend_handles_labels()
    h2, l2 = axs2.get_legend_handles_labels()
    axs.legend(h1+h2, l1+l2, loc = 'lower right')
    axs.grid('major')
    axs.set_xlabel('z ($\\mu$m)')
    axs.set_ylabel('U/$k_{\\rm B}$ ($\\mu$K)')
    axs2.set_ylabel('e$^{-U/k_{\\rm B}T}$')
    axs.set_title('T = {:.1f}'.format(T*1e6)+' $\\mu$K, N$_{\\rm grid}$ = '+' {:.0f}'.format(N))
    
    axs = fig.add_subplot(gs[1,0]) 
    idx_z = np.abs(Z[0,0,:]-z0).argmin()
    idx_y = np.abs(Y[0,:,0]-0).argmin()
    axs.plot(X[:, idx_y,idx_z]*1e6, Uval[:, idx_y,idx_z]*1e6/thermo.kB, '-', label = 'U at (x,0,{:.0f})'.format(z0_rd*1e6), \
             color = clrPts[0], linewidth = 2)
    axs.set(ylim = (plot_min,plot_max))
    axs2 = axs.twinx()
    axs2.plot(X[:, idx_y,idx_z]*1e6, integrand[:, idx_y,idx_z], '-', label = 'e$^{-U/k_{\\rm B}T}$'+' at (x,0,{:.0f})'.format(z0_rd*1e6),\
              color = clrPts[1], linewidth = 2)
    axs2.set(ylim = (0,1))
    h1, l1 = axs.get_legend_handles_labels()
    h2, l2 = axs2.get_legend_handles_labels()
    axs.legend(h1+h2, l1+l2, loc = 'upper right')
    axs.grid('major')
    axs.set_xlabel('x ($\\mu$m)')
    axs.set_ylabel('U/$k_{\\rm B}$ ($\\mu$K)')
    axs2.set_ylabel('e$^{-U/k_{\\rm B}T}$')
    
    axs = fig.add_subplot(gs[2,0]) 
    axs.plot(l*1e6,V*1e6, '.', label = 'V', color = clrPts[0], markersize = 10)
    axs2 = axs.twinx()
    axs2.plot(l*1e6,S/thermo.kB,'.', label = 'S/N', color = clrPts[1], markersize = 10)
    h1, l1 = axs.get_legend_handles_labels()
    h2, l2 = axs2.get_legend_handles_labels()
    axs.legend(h1+h2, l1+l2, loc = 'lower right')
    axs.grid('major')
    axs.set_xlabel('L ($\\mu$m)')
    axs.set_ylabel('V (cm$^3$)')
    axs2.set_ylabel('S/N/k$_{\\rm B}$ ($\\mu$K)')
    
    return fig
    
# %% Quick Inputs #############################################################
Tdisp         = 50e-6      # Temperature at which convergence is to be evaluated [K]
N_grid        = 201        # Spatial Grid Size [odd int]
size_grid     = 30          # Spatial Grid extent [fraction of IR beam waist w0_rd]
# %% Other Inputs #############################################################
N_Tscale      = 100        # Temperature Scale grid size [int]
Tscale_start  = 0.8e-6     # Start of the temperature scale [K]
Tscale_end    = 100e-6     # End of the temperature scale [K]
# Quadrupole Magnetic Trap Properties =========================================
Bdash         = 160        # Field Gradient [G/cm]
# IR Dipole Trap Properties ===================================================
P0_rd         = 50e-6      # Total Power [W]
trapFrqYZ0_rd = 20         # Trap Frequency at P0_rd [Hz]
odtWL_rd      = 1064e-9    # Wavelength [m]
w0_rd         = 65*1e-6    # Beam waist size [m]
P_rd          = 0.0035     # Power [W]
z0_rd         = -1*w0_rd   # Beam waist vertical position [m]

# %% Initialize ###############################################################    
trapFrqYZ_rd = odt.freqVSpower(trapFrqYZ0_rd, P0_rd, P_rd)
U0  = odt.odtU_depth(na.m, nu = [trapFrqYZ_rd, trapFrqYZ_rd*odtWL_rd/np.pi/w0_rd], \
        r = [w0_rd, np.pi*w0_rd**2/odtWL_rd])
def U_total(bdash, X,Y,Z):
    U_rd = odt.odtU_redDet_singleBeam(U0, w0_rd, odtWL_rd, 0,0,z0_rd, X, Y, Z)
    U   = U_rd \
          + odt.odtU_quadMagnetic(bdash, na.fieldSplit, X,Y,Z) \
          + odt.gravityU(Z, na.m)
    U = U - np.min(U)
    return U
x,y,z = gridMnp._uniform1Dgrids(size_grid*w0_rd, N_grid)           

# %% Evaluate Convergence #####################################################
l_max      = np.linspace(2*np.mean(np.diff(x)),size_grid*w0_rd,int((N_grid-1)/2))
VatTdisp   = np.zeros(np.shape(l_max))
SatTdisp   = np.zeros(np.shape(l_max))

for ii in range(len(l_max)):
    X,Y,Z, x_temp, y_temp, z_temp = gridMnp._make3Dgrid(x, y, z, l_max[ii], l_max[ii], l_max[ii])
    U   = U_total(Bdash, X, Y, Z)
    odtPropScale = thermo.thermoODT( U, [x_temp,y_temp,z_temp], na.m, \
             T_start = Tscale_start, T_end = Tscale_end , N_T = N_Tscale)  
    odtPropScale.entropyPerAtom()
    Tdisp, VatTdisp[ii], _, _, SatTdisp[ii] = odtPropScale.propsAtT0(Tdisp)

# %% Evaluate potential Energy landscape for plot #############################
X,Y,Z, _,_,_ = gridMnp._make3Dgrid(x, y, z, size_grid*w0_rd, size_grid*w0_rd, size_grid*w0_rd)
Uval = U_total(Bdash, X, Y, Z)
integrand = np.exp(-Uval/thermo.kB/Tdisp) 

# %% Plot #####################################################################
fig = plotConv(1,X,Y,Z,Uval,Tdisp,N_grid,z0_rd,l_max,SatTdisp,VatTdisp)
fig.savefig(paths.savDir+'convgTest_{:.0f}uK_{:.0f}.png'.format(Tdisp*1e6,N_grid))
