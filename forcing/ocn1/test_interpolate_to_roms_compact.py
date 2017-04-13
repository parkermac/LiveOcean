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

This version of the code is meant to prototype the functions we will use
in make_forcing_main.py.

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


import zfun
import zrfun
import pickle
import numpy as np
import time

import Ofun

from importlib import reload
reload(Ofun)
reload(zfun)

in_dir = Ldir['LOogf_fd']

# get grid and S info
G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
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

def get_xyr(G, vn):
    # ROMS lon, lat arrays
    if vn in ['ssh', 'theta', 's3d']:
        xr = G['lon_rho']
        yr = G['lat_rho']
    elif vn in ['ubar', 'u3d']:
        xr = G['lon_u']
        yr = G['lat_u']
    elif vn in ['vbar', 'v3d']:
        xr = G['lon_v']
        yr = G['lat_v']
    else:
        print('Unknown variable name for get_xyr: ' + vn)
    return xr, yr

def get_zr(G, S, vn):
    # get z on the ROMS grids
    h = G['h']
    if vn in ['theta', 's3d']:
        zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    elif vn in ['u3d']:    
        xru, yru = get_xyr(G, 'ubar')
        hu = zfun.interp_scattered_on_plaid(G['lon_u'], G['lat_u'],
                    G['lon_rho'][0,:], G['lat_rho'][:,0], h, exnan=False)
        hu = np.reshape(hu, G['lon_u'].shape)
        zr = zrfun.get_z(hu, 0*hu, S, only_rho=True)    
    elif vn in ['v3d']:    
        hv = zfun.interp_scattered_on_plaid(G['lon_v'], G['lat_v'],
                    G['lon_rho'][0,:], G['lat_rho'][:,0], h, exnan=False)
        hv = np.reshape(hv, G['lon_v'].shape)
        zr = zrfun.get_z(hv, 0*hv, S, only_rho=True)
    else:
        print('Unknown variable name for get_zr: ' + vn)
    return zr

def get_interpolated(G, S, b, lon, lat, z, N):
    # start input dict
    c = dict()
    # 2D fields
    for vn in ['ssh', 'ubar', 'vbar']:
        xr, yr = get_xyr(G, vn)
        Mr, Lr = xr.shape
        u = b[vn]
        uu = zfun.interp_scattered_on_plaid(xr, yr, lon, lat, u)
        uu = np.reshape(uu, xr.shape)
        if np.isnan(uu).any():
            print('Warning: nans in output array for ' + vn)
        else:
            c[vn] = uu
    # 3D fields
    for vn in ['theta', 's3d', 'u3d', 'v3d']:
        xr, yr = get_xyr(G, vn)
        zr = get_zr(G, S, vn)
        Nr, Mr, Lr = zr.shape
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
            print('Warning: nans in output array for ' + vn)
        else:
            c[vn] = vv
    return c

tt0 = time.time()
c = get_interpolated(G, S, b, lon, lat, z, N)   
print('Interpolation to ROMS took %0.1f seconds' % (time.time() - tt0))





