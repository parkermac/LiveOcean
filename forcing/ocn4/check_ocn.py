# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 07:14:06 2016

@author: PM5

Code to look at the ocn files on a given day.

Only set up to work on mac.

"""

gridname='cas6'
tag='v3'
date_string = '2020.08.12'

vn = 'salt' # roms variable name to plot
cmap = 'rainbow' # default colormap


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

lonr = dsg['lon_rho'][:]
latr = dsg['lat_rho'][:]
maskr = dsg['mask_rho'][:]

lonu = dsg['lon_u'][:]
latu = dsg['lat_u'][:]
masku = dsg['mask_u'][:]

lonv = dsg['lon_v'][:]
latv = dsg['lat_v'][:]
maskv = dsg['mask_v'][:]

dsg.close()

indir = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/ocn4/'
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

#plt.close()


def set_box(ax):
    ax.set_xlim(-130.5, -122)
    ax.set_ylim(41.5,52.5)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18,7), squeeze=False)

if vn == 'salt':
    hvn = 's3d'
    vmin = 24
    vmax = 34
elif vn == 'temp':
    hvn = 't3d'
    vmin = 5
    vmax = 12
elif vn == 'zeta':
    hvn = 'ssh'
    vmin = 0
    vmax = .4
elif vn == 'ubar':
    hvn = 'ubar'
    vmin = -.2
    vmax = .2
elif vn == 'vbar':
    hvn = 'vbar'
    vmin = -.2
    vmax = .2

# extract variables to plot
if vn in ['zeta', 'ubar', 'vbar']:
    # 2D
    vxfh = xfh[hvn][:,:]
    try:
        vfh = fh[hvn][:,:]
    except KeyError: # no vbar in fh; it is created during extrapolation
        vfh = 0 * vxfh
    vroms = ds[vn][0, :, :]
else:
    # 3D
    vfh = fh[hvn][-1,:,:]
    vxfh = xfh[hvn][-1,:,:]
    vroms = ds[vn][0, -1, :, :]
    
ax = axes[0,0]
cs = ax.pcolormesh(lon, lat, vfh, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.dar(ax)
pfun.add_coast(ax)
set_box(ax)
ax.set_title('HYCOM original: ' + hvn)

ax = axes[0,1]
cs = ax.pcolormesh(lon, lat, vxfh, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.dar(ax)
pfun.add_coast(ax)
set_box(ax)
ax.set_title('HYCOM extrapolated: ' + hvn)

ax = axes[0,2]

if vn == 'ubar':
    Lon = lonu; Lat = latu; Mask = masku
elif vn == 'vbar':
    Lon = lonv; Lat = latv; Mask = maskv
else:
    Lon = lonr; Lat = latr; Mask = maskr
vroms[Mask==0] = np.nan

cs = ax.pcolormesh(Lon, Lat, vroms, cmap=cmap, vmin=vmin, vmax=vmax)
fig.colorbar(cs)
pfun.dar(ax)
pfun.add_coast(ax)
ax.text(.95, .05, date_string,
    horizontalalignment='right', transform=ax.transAxes,
    fontweight='bold')
set_box(ax)
ax.set_title('ROMS ' + Ldir['gtag'] + ': ' + vn)

plt.show()

