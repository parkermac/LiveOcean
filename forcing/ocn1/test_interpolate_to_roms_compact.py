#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 16:59:30 2017

@author: PM5

This is a test of taking a set of HyCOM fields and interpolating them
to a ROMS grid.  We do it in two steps, first interpolating horizontally
to the ROMS lon, lat, and then interpolating vertically to the ROMS z.

By making maximum use of flattening and reshaping to avoid loops the
performance is excellent (and gives results that are identical to the versions
that use loops).  We achieve about a 25x speed up.
"""

import os
import sys
fpth = os.path.abspath('../../alpha')
if fpth not in sys.path:
    sys.path.append(fpth)
import Lfun
Ldir = Lfun.Lstart('cas1', 'base')
Ldir['date_string'] = '2013.01.01'
Ldir['LOogf_fd'] = (Ldir['LOo'] + Ldir['gtag'] + '/f' + Ldir['date_string']
        + '/ocn1/Data/')

import matplotlib.pyplot as plt

import zfun
import zrfun
import pickle
import numpy as np
import time

import Ofun
from importlib import reload
reload(Ofun)

in_dir = Ldir['LOogf_fd']

# get grid and S info
[G] = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', getS=False, getT=False)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

a = os.listdir(in_dir)
aa = [item for item in a if item[:3]=='xfh']


#%% do the extrapolated using function calls

# HyCOM grid info
lon, lat, z, L, M, N, X, Y = Ofun.get_coords(in_dir)

# load a dict of hycom fields
in_fn = in_dir + aa[0]
b = pickle.load(open(in_fn, 'rb'))

# ROMS lon, lat arrays
xr = G['lon_rho']
yr = G['lat_rho']
# ROMS grid sizes
Nr = S['N']
Mr, Lr = xr.shape

c = dict()

# 2D fields
for vn in ['ssh', 'ubar', 'vbar']:
    u = b[vn]
    uu = zfun.interp_scattered_on_plaid(xr, yr, lon, lat, u)
    uu = np.reshape(uu, xr.shape)
    c[vn] = uu

# choose a 3D variable to work on
for vn in ['theta', 's3d', 'u3d', 'v3d']:
    v = b[vn]
    # create an intermediate array, which is on the ROMS lon, lat grid
    # but has the HyCOM vertical grid
    # get interpolants
    xi0, xi1, xf = zfun.get_interpolant(xr, lon, extrap_nan=True)
    yi0, yi1, yf = zfun.get_interpolant(yr, lat, extrap_nan=True)
    # bi linear interpolation
    u00 = v[:,yi0,xi0]
    u10 = v[:,yi1,xi0]
    u01 = v[:,yi0,xi1]
    u11 = v[:,yi1,xi1]
    vi = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)
    vi = vi.reshape((N, Mr, Lr))
    # get z on the ROMS grid    
    h = G['h']
    zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    # make interpolants to go from the HyCOM vertical grid to the
    # ROMS vertical grid
    I0, I1, FR = zfun.get_interpolant(zr, z, extrap_nan=True)
    zrs = zr.shape
    
    vif = vi.flatten()
    LL = np.tile(np.arange(Lr), Mr*Nr)
    MM = np.tile(np.repeat(np.arange(Mr), Lr), Nr)
    f0 = LL + MM*Lr + I0*Mr*Lr
    f1 = LL + MM*Lr + I1*Mr*Lr
    vv = (1-FR)*vif[f0] + FR*vif[f1]
    vv = vv.reshape(zrs)

if np.isnan(vv).any():
    print('Warning: nans in output array')
    
#%% plotting

plt.close('all')

# 2D fields
fig, axes = plt.subplots(2, 3, figsize=(13,8))

cmap = 'rainbow'

ax = axes[0, 0]
pc = ax.pcolormesh(lon, lat, v[-1,:,:], cmap=cmap)
fig.colorbar(pc, ax=ax)
A = ax.axis()

ax = axes[0, 1]
pc = ax.pcolormesh(xr, yr, vi[-1,:,:], cmap=cmap)
fig.colorbar(pc, ax=ax)
ax.axis(A)

ax = axes[0, 2]
pc = ax.pcolormesh(xr, yr, vv[-1,:,:], cmap=cmap)
fig.colorbar(pc, ax=ax)
ax.axis(A)

ax = axes[1, 0]
pc = ax.pcolormesh(lon, z, v[:,int(M/2),:], cmap=cmap)
fig.colorbar(pc, ax=ax)
A = ax.axis()

ax = axes[1, 1]
pc = ax.pcolormesh(xr[int(Mr/2),:], z, vi[:,int(Mr/2),:], cmap=cmap)
fig.colorbar(pc, ax=ax)
ax.axis(A)

ax = axes[1, 2]
pc = ax.pcolormesh(xr[int(Mr/2),:], zr[:,int(Mr/2),:], vv[:,int(Mr/2),:], cmap=cmap)
fig.colorbar(pc, ax=ax)
ax.axis(A)

plt.show()




