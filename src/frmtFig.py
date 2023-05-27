"""
###############################################################################
Code for formatting figures
###############################################################################
Created:    Swarnav Banik  on  May 12, 2023 
"""

def frmtFig(mpl, plt, FS_title = 20, FS_tickLabel = 20, FS_axisLabel = 20):
    # Colors ##################################################################
    clrBckg = (8/100, 11/100, 17.6/100)
    clrText = (1,1,1)
    clrPts1 = (90/100, 47/100, 27/100)
    clrPts2 = (0/100, 50/100, 99/100)
    clrPts3 = (74/100, 76/100, 75/100)
    clrPts4 = (31/100, 74/100, 55/100)
    clrPts  = [clrPts1, clrPts2, clrPts3, clrPts4]

    mpl.rcdefaults()
    # Set figure font sizes ###################################################
    mpl.rcParams['axes.labelsize']  = FS_axisLabel
    mpl.rcParams['xtick.labelsize'] = FS_tickLabel
    mpl.rcParams['ytick.labelsize'] = FS_tickLabel
    mpl.rcParams['legend.fontsize'] = FS_axisLabel
    mpl.rcParams['axes.titlesize']  = FS_title
    plt.rc('figure', titlesize = FS_title)
    plt.rcParams['text.color'] = clrText
    # Set the figure color ####################################################
    mpl.rcParams['axes.facecolor'] = clrBckg
    mpl.rcParams['axes.edgecolor'] = clrText
    mpl.rcParams['axes.labelcolor'] = clrText
    mpl.rcParams['xtick.color'] = clrText
    mpl.rcParams['ytick.color'] = clrText
    mpl.rcParams['figure.facecolor'] = clrBckg
    # Set the grid ############################################################
    mpl.rcParams['grid.color'] = clrText
    mpl.rcParams['grid.alpha'] = 0.4
    mpl.rcParams['grid.linewidth'] = 0.8
    # Set Legend properties ###################################################
    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['legend.title_fontsize'] = FS_axisLabel
        
    mpl.rcParams.update(mpl.rcParams)
    
    return clrPts, mpl, plt