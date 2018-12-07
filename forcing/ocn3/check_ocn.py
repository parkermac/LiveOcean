# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 07:14:06 2016

@author: PM5

Code to look at the ocn files on a given day.

Only set up to work on mac.

"""

gridname='cas4'
tag='v2'
date_string = '2018.12.07'

# setup

import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname=gridname, tag=tag)
import zfun
import zrfun

pth = os.path.abspath('../../plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle

# get files

grid_fn = Ldir['grid'] + 'grid.nc'
dsg = nc.Dataset(grid_fn)
lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]
maskr = dsg['mask_rho'][:]
dsg.close()

indir = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/ocn3/'
in_fn = (indir + 'ocean_ini.nc')
try:
    ds = nc.Dataset(in_fn)
except OSError:
    pass
    
in_fn_coords = (indir + 'Data/coord_dict.p')
in_fn_fh = (indir + 'Data/fh' + date_string + '.p')
in_fn_xfh = (indir + 'Data/xfh' + date_string + '.p')
#
coord_dict = pickle.load(open(in_fn_coords, 'rb'))
lon = coord_dict['lon']
lat = coord_dict['lat']
#
fh = pickle.load(open(in_fn_fh, 'rb'))
xfh = pickle.load(open(in_fn_xfh, 'rb'))

#%% plotting

plt.close()


def set_box(ax):
    ax.set_xlim(-127.5, -122)
    ax.set_ylim(42.5,50.5)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(13,7), squeeze=False)

if False:
    hvn = 's3d'
    vn = 'salt'
    vmin = 15
    vmax = 34
else:
    hvn = 't3d'
    vn = 'temp'
    vmin = 5
    vmax = 12

ax = axes[0,0]
v = fh[hvn][-1,:,:]
cs = ax.pcolormesh(lon, lat, v, cmap='rainbow', vmin=vmin, vmax=vmax)
pfun.dar(ax)
pfun.add_coast(ax)
set_box(ax)
ax.set_title('HYCOM original: ' + hvn)

ax = axes[0,1]
v = xfh[hvn][-1,:,:]
cs = ax.pcolormesh(lon, lat, v, cmap='rainbow', vmin=vmin, vmax=vmax)
pfun.dar(ax)
pfun.add_coast(ax)
set_box(ax)
ax.set_title('HYCOM extrapolated: ' + hvn)

try:
    ax = axes[0,2]
    v = ds[vn][0, -1, :, :]
    v[maskr==0] = np.nan
    cs = ax.pcolormesh(lonp, latp, v[1:-1,1:-1], cmap='rainbow', vmin=vmin, vmax=vmax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.text(.95, .05, date_string,
        horizontalalignment='right', transform=ax.transAxes,
        fontweight='bold')
    set_box(ax)
    ax.set_title('ROMS ' + Ldir['gtag'] + ': ' + vn)
except:
    pass

plt.show()

