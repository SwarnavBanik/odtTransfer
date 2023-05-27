from src.frmtFig import frmtFig 
import matplotlib as mpl
import matplotlib.pyplot as plt
clrPts, mpl, plt = frmtFig(mpl, plt, FS_title = 18, FS_tickLabel = 18, FS_axisLabel = 18)

import numpy as np
import src.rb as rb
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
              label = '(x = {:.0f}, y = {:.0f})'.format(X[idx_x,0,0]*1e6,Y[0,idx_y,0]*1e6),\
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
              label = '(y = {:.0f}, z = {:.0f})'.format(Y[0,idx_y,0]*1e6,Z[0,0,idx_z]*1e6),\
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
                 ', y,z) /k$_{\\rm B}$ ($\\mu$K)', orientation='horizontal')
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
                 ') /k$_{\\rm B}$ ($\\mu$K)', orientation='horizontal')
    axs.set_ylabel('x ($\\mu$m)')
    axs.set_xlabel('y ($\\mu$m)')
    axs.set(xlim = (-plot_l, +plot_l))
    axs.set(ylim = (-plot_l, +plot_l))
    
    return fig

def plotThermo(odt_0, odt_f, T,V,eta,S, figNo = 1):
    
    D = 1/eta
    fig = plt.figure(figNo, figsize = (12*1,3*6))
    gs = GridSpec(3,1)
    
    axs1 = fig.add_subplot(gs[0,0]) 
    axs1.loglog(odt_0.T*1e6, odt_0.V*1e6,'.', label = 'U$_{\\rm i}$', linewidth = 5,\
               color = clrPts[0], markersize = 10)
    axs1.loglog(odt_f.T*1e6, odt_f.V*1e6,'.',  label = 'U$_{\\rm f}$', linewidth = 5,\
               color = clrPts[1], markersize = 10)
    axs1.loglog(T*1e6, V*1e6,'--', label = 'Isentropic Process', \
                linewidth = 3, color = clrPts[3], markersize = 10)
    axs1.set_ylabel('V (cm$^3$)', fontsize = 22)
    h1, l1 = axs1.get_legend_handles_labels()
    axs1.legend(h1, l1, loc = 'lower right')
    axs1.grid('major')
    
    axs2 = fig.add_subplot(gs[1,0]) 
    axs2.semilogx(odt_0.T*1e6, odt_0.S/thermo.kB,'.', label = 'U$_{\\rm i}$', \
                 linewidth = 5, color = clrPts[0], markersize = 10)
    axs2.semilogx(odt_f.T*1e6, odt_f.S/thermo.kB,'.', label = 'U$_{\\rm f}$', \
                 linewidth = 5, color = clrPts[1], markersize = 10)
    axs2.semilogx(T*1e6, S/thermo.kB,'--', label = 'Isentropic Process', \
                  linewidth = 3, color = clrPts[3], markersize = 10)    
    axs2.set_ylabel('S/Nk$_{\\rm B}$ ($\\mu$K)', fontsize = 22)
    axs2.grid('major')
    h1, l1 = axs2.get_legend_handles_labels()
    axs2.legend(h1, l1, loc = 'upper left')
        
    axs3 = fig.add_subplot(gs[2,0])    
    axs3.loglog(odt_0.T*1e6, 1/odt_0.eta,'.', label = 'U$_{\\rm i}$', \
                 linewidth = 5, color = clrPts[0], markersize = 10)
    axs3.loglog(odt_f.T*1e6, 1/odt_f.eta,'.', label = 'U$_{\\rm f}$', \
                 linewidth = 5, color = clrPts[1], markersize = 10)
    axs3.loglog(T*1e6, D,'--', label = 'Isentropic Process', \
                  linewidth = 3, color = clrPts[3], markersize = 10)
    axs3.set_xlabel('T ($\\mu$K)', fontsize = 22)
    axs3.set_ylabel('D/N ($\\mu$K)', fontsize = 22)
    h1, l1 = axs3.get_legend_handles_labels()
    axs3.legend(h1, l1, loc = 'lower left')
    axs3.grid('major')    
    return fig

# %% Quick Inputs #############################################################
N_is          = 50         # No. of steps for isentropic process [int]
T0            = 50e-6      # Initial temperature [K]
# %% Other Inputs #############################################################
N_grid        = 201        # Spatial Grid Size [odd int]
size_grid     = 30         # Spatial Grid extent
N_Tscale      = 100        # Temperature Scale grid size [int]
Tscale_start  = 0.8e-6     # Start of the temperature scale [K]
Tscale_end    = 100e-6     # End of the temperature scale [K]
# Quadrupole Magnetic Trap Properties =========================================
Bdash_0       = 160        # Initial Field Gradient [G/cm]
Bdash_f       = 31         # Firbl Field Gradient [G/cm]
# IR Dipole Trap Properties ===================================================
P0_rd         = 50e-6      # Total Power [W]
trapFrqYZ0_rd = 20         # Trap Frequency at P0_rd [Hz]
odtWL_rd      = 1064e-9    # Wavelength [m]
w0_rd         = 65*1e-6    # Beam waist size [m]
P_rd          = 0.0035     # Power [W]
z0_rd         = -1*w0_rd   # Beam waist vertical position [m]

# %% Initialize ###############################################################    
Bdash = np.linspace(Bdash_0, Bdash_f, N_is)              # Gradient [G/cm]
# Spatial grid ================================================================
x,y,z = gridMnp._uniform1Dgrids(size_grid*w0_rd, N_grid)   
X,Y,Z, x,y,z = gridMnp._make3Dgrid(x, y, z, size_grid*w0_rd, size_grid*w0_rd, size_grid*w0_rd)
# Trap Potential ==============================================================
trapFrqYZ_rd = odt.freqVSpower(trapFrqYZ0_rd, P0_rd, P_rd)
U0  = odt.odtU_depth(rb.m, nu = [trapFrqYZ_rd, trapFrqYZ_rd*odtWL_rd/np.pi/w0_rd], \
        r = [w0_rd, np.pi*w0_rd**2/odtWL_rd])
U_rd = odt.odtU_redDet_singleBeam(U0, w0_rd, odtWL_rd, 0,0,z0_rd, X, Y, Z)
def U_total(bdash):
    U   = U_rd \
          + odt.odtU_quadMagnetic(bdash, rb.fieldSplit, X,Y,Z) \
          + odt.gravityU(Z, rb.m)
    U = U - np.min(U)
    return U
# Thermodyrbmic properties v/s Tscale =========================================
def odtPropsvsTscale_atBdash(bdash, figNo = 1, plot = False):
    if np.size(bdash) != 1:
        raise Exception('adibaticExpn:: SvsTscale_atBdash: bdash needs to be of size 1')
    U   = U_total(bdash)
    if plot:
        fig = plot3DU(U,X,Y,Z, z0 = z0_rd, figNo = figNo, plot_min= -20, plot_max = 150, plot_l = [])
        
    odtPropScale = thermo.thermoODT( U, [x,y,z], rb.m, \
             T_start = Tscale_start, T_end = Tscale_end , N_T = N_Tscale) 
    odtPropScale.entropyPerAtom()
    if plot:
        return odtPropScale, fig
    else:
        return odtPropScale

# %% Evaluate the initial and Firbl thermodynmic properties ###################
odt_0, fig1 = odtPropsvsTscale_atBdash(Bdash_0, figNo = 1, plot = True)
odt_f, fig2 = odtPropsvsTscale_atBdash(Bdash_f, figNo = 2, plot = True)
fig1.savefig(paths.savDir+'adiabaticExpn_U_initial.png')
fig2.savefig(paths.savDir+'adiabaticExpn_U_firbl.png')

# %% Evaluate the initial conditions ##########################################
T0, V0, _, _, S0 = odt_0.propsAtT0(T0)
print('The trap starts from T = {:.1f} uK and S/N = {:.1f} kB'.format(T0*1e6, S0/thermo.kB))

# %% Evaluate the isentropic process ##########################################
T    = np.zeros(np.shape(Bdash))
S    = np.zeros(np.shape(Bdash))
V    = np.zeros(np.shape(Bdash))
eta  = np.zeros(np.shape(Bdash))
for ii in range(len(Bdash)):
    odt_intr = odtPropsvsTscale_atBdash(Bdash[ii], plot = False)
    T[ii], V[ii], eta[ii], _, S[ii] = odt_intr.propsAtS0(S0)
# %% Plot the thermodyrbmic properties ########################################
fig = plotThermo(odt_0, odt_f, T,V,eta,S, figNo = 3)
fig.savefig(paths.savDir+'adiabaticExpn_{:.0f}uK_{:.0f}.png'.format(T0*1e6,N_is))

    


