"""
###############################################################################
Code for manipulating 3D grid
###############################################################################
Created:    Swarnav Banik  on  May 26, 2023 
"""

import numpy as np

def _uniform1Dgrid(l_max, N):
    if l_max <= 0:
        raise Exception('gridMnp:: _uniform1Dgrid: l_max should be positive.')
    return np.linspace(-l_max, +l_max, N)

def _uniform1Dgrids(l_max, N):
    if l_max <= 0:
        raise Exception('gridMnp:: _uniform1Dgrid: l_max should be positive.')
    r = _uniform1Dgrid(l_max, N)
    return r,r,r
    
def _truncateGrid(A,l_max):
    if np.ndim(A) != 1:
        raise Exception('gridMnp:: _truncateGrid: A should be an array of dimension 1.')
    idx_start = np.abs(A + l_max).argmin()
    idx_end   = np.abs(A - l_max).argmin()
    return A[idx_start:idx_end+1]
    
def _make3Dgrid(x, y, z, x_max, y_max, z_max):
    if np.ndim(x) != 1 or np.ndim(y) != 1 or np.ndim(z) != 1:
        raise Exception('convergenceTest:: _make3Dgrid: x, y, z should be an arrays of dimension 1.')
    if np.size(x_max) != 1 or np.size(y_max) != 1 or np.size(z_max) != 1:
        raise Exception('convergenceTest:: _make3Dgrid: x_max, y_max, and z_max should be of size 1.')
    x = _truncateGrid(x,x_max)
    y = _truncateGrid(y,y_max)
    z = _truncateGrid(z,z_max)
    X, Y, Z = np.meshgrid(x, y, z, indexing = 'ij') 
    return X,Y,Z, x,y,z
