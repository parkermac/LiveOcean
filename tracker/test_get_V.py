# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:23:53 2016

@author: PM5

Code to test the routine for getting field values at arbitrary points.

"""

#%% setup
import numpy as np
import netCDF4 as nc4

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
import zfun
reload(zfun)

#%%

fn = ('/Users/PM5/Documents/LiveOcean_roms/output/cascadia1_base_lobio1/'
    + 'f2015.08.15/ocean_his_0002.nc')

# get basic info
G, S = zfun.get_basic_info(fn, getT=False)

# make vectors to feed to interpolant maker
R = dict()
R['rlonr'] = G['lon_rho'][0,:].squeeze()
R['rlatr'] = G['lat_rho'][:,0].squeeze()
R['rlonu'] = G['lon_u'][0,:].squeeze()
R['rlatu'] = G['lat_u'][:,0].squeeze()
R['rlonv'] = G['lon_v'][0,:].squeeze()
R['rlatv'] = G['lat_v'][:,0].squeeze()
R['rcsr'] = S['Cs_r'][:]
R['rcsw'] = S['Cs_w'][:]

# these lists are used internally to get other variables as needed
vn_list = ['salt', 'temp', 'zeta', 'h', 'u', 'v', 'w',
                 'Uwind', 'Vwind']

ds = nc4.Dataset(fn)

plon = np.linspace(-124.5, - 123.5, 10)
plat = 46.5 * np.ones_like(plon)
pcs = -0.5 * np.ones_like(plon)

# THE FUNCTION
# would be defined as:
# def get_V(vn_list, ds, plon, plat, pcs, R):

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# get interpolant arrays
i0lon_d = dict()
i1lon_d = dict()
frlon_d = dict()
i0lat_d = dict()
i1lat_d = dict()
frlat_d = dict()
for gg in ['r', 'u', 'v']:
    exn = False
    i0lon_d[gg], i1lon_d[gg], frlon_d[gg] = zfun.get_interpolant(
            plon, R['rlon'+gg], extrap_nan=exn)
    i0lat_d[gg], i1lat_d[gg], frlat_d[gg] = zfun.get_interpolant(
            plat, R['rlat'+gg], extrap_nan=exn)
i0csr, i1csr, frcsr = zfun.get_interpolant(pcs, R['rcsr'], extrap_nan=exn)
i0csw, i1csw, frcsw = zfun.get_interpolant(pcs, R['rcsw'], extrap_nan=exn)
NV = len(vn_list)
NP = len(plon)
# get interpolated values
V = np.nan * np.ones((NP,NV))
vcount = 0
for vn in vn_list:
    if vn in ['w']:
        i0cs = i0csw
        i1cs = i1csw
        frcs = frcsw
    else:
        i0cs = i0csr
        i1cs = i1csr
        frcs = frcsr
    if vn in ['salt','temp','zeta','h','Uwind','Vwind', 'w']:
        gg = 'r'
    elif vn in ['u']:
        gg = 'u'
    elif vn in ['v']:
        gg = 'v'
    i0lat = i0lat_d[gg]
    i1lat = i1lat_d[gg]
    frlat = frlat_d[gg]
    i0lon = i0lon_d[gg]
    i1lon = i1lon_d[gg]
    frlon = frlon_d[gg]
    # get the data field and put nan's in masked points
    v0 = ds.variables[vn][:].squeeze()
    try:
        vv = v0.data
        vv[v0.mask] = np.nan
    except AttributeError:
        # it is not a masked array
        vv = v0
    if vn in ['salt','temp','u','v','w']:
        # Get just the values around our particle positions.
        # each row in VV corresponds to a "box" around a point
        VV = np.nan* np.ones((NP, 8))
        VV[:,0] = vv[i0cs, i0lat, i0lon]
        VV[:,1] = vv[i0cs, i0lat, i1lon]
        VV[:,2] = vv[i0cs, i1lat, i0lon]
        VV[:,3] = vv[i0cs, i1lat, i1lon]
        VV[:,4] = vv[i1cs, i0lat, i0lon]
        VV[:,5] = vv[i1cs, i0lat, i1lon]
        VV[:,6] = vv[i1cs, i1lat, i0lon]
        VV[:,7] = vv[i1cs, i1lat, i1lon]
        # Work on edge values.  If all in a box are masked
        # then that row will be nan's, and also:
        if vn in ['u', 'v', 'w']:
            # set all velocities to zero if any in the box are masked
            VV[np.isnan(VV).any(axis=1), :] = 0
        elif vn in ['salt','temp']:
            # set all tracers to their average if any in the box are masked
            newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,8))
            mask = np.isnan(VV)
            VV[mask] = newval[mask]
        # now do the interpolation in each box
        vl = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
            + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
        vu = ( (1-frlat)*((1-frlon)*VV[:,4] + frlon*VV[:,5])
            + frlat*((1-frlon)*VV[:,6] + frlon*VV[:,7]) )
        v = (1-frcs)*vl + frcs*vu
    elif vn in ['zeta','Uwind','Vwind', 'h']:
        VV = np.nan* np.ones((NP, 4))
        VV[:,0] = vv[i0lat, i0lon]
        VV[:,1] = vv[i0lat, i1lon]
        VV[:,2] = vv[i1lat, i0lon]
        VV[:,3] = vv[i1lat, i1lon]
        newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
        mask = np.isnan(VV)
        VV[mask] = newval[mask]
        v = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
            + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
    V[:,vcount] = v
    vcount += 1

