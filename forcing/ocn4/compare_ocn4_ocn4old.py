# -*- coding: utf-8 -*-
"""
Code to compare the results of ocn4 to ocn4old, to make sure
we are doing the interpolation well enough.

"""

gridname='cas6'
tag='v3'
date_string = '2017.04.20'

layer = 'surface' # 'surface' or 'bottom'

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
ds = nc.Dataset(in_fn)
    
indir = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/ocn4old/'
in_fn = (indir + 'ocean_ini.nc')
dso = nc.Dataset(in_fn)

vlim_dict = {'salt': (15,34,1),
            'temp': (5,12,1),
            'u': (-.4,.4,.04),
            'v': (-.4,.4,.04),
            'zeta': (0,.4,.02),
            'ubar': (-.2,.2,.02),
            'vbar': (-.2,.2,.02)}

#%% plotting
aa = [-130, -122, 42, 52]
plt.close('all')

            
for vn in vlim_dict.keys():
    
    vmin = vlim_dict[vn][0]
    vmax = vlim_dict[vn][1]
    dv = vlim_dict[vn][2]

    # extract variables to plot
    if vn in ['zeta', 'ubar', 'vbar']:
        # 2D
        vv = ds[vn][0, :, :]
        vvo = dso[vn][0, :, :]
    else:
        # 3D
        if layer == 'surface':
            vv = ds[vn][0, -1, :, :]
            vvo = dso[vn][0, -1, :, :]
        elif layer == 'bottom':
            vv = ds[vn][0, 0, :, :]
            vvo = dso[vn][0, 0, :, :]
    
    if vn in ['u', 'ubar']:
        Lon = lonu; Lat = latu; Mask = masku
    elif vn in ['v', 'vbar']:
        Lon = lonv; Lat = latv; Mask = maskv
    else:
        Lon = lonr; Lat = latr; Mask = maskr
    vv[Mask==0] = np.nan
    vvo[Mask==0] = np.nan
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,7), squeeze=False)

    ax = axes[0,0]
    cs = ax.pcolormesh(Lon, Lat, vvo, cmap='rainbow', vmin=vmin, vmax=vmax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title(layer.title() + ' ocn4old')
    ax.text(.97, .15, vn,
        horizontalalignment='right', transform=ax.transAxes,
        fontsize=18, fontweight='bold')
    ax.text(.97, .10, Ldir['gtag'],
        horizontalalignment='right', transform=ax.transAxes)
    ax.text(.97, .05, date_string,
        horizontalalignment='right', transform=ax.transAxes)
    ax.axis(aa)

    ax = axes[0,1]
    cs = ax.pcolormesh(Lon, Lat, vv, cmap='rainbow', vmin=vmin, vmax=vmax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('ocn4')
    ax.axis(aa)

    ax = axes[0,2]
    cs = ax.pcolormesh(Lon, Lat, vv - vvo, cmap='bwr', vmin=-dv, vmax=dv)
    fig.colorbar(cs)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('ocn4 - ocn4old')
    ax.axis(aa)

plt.show()

