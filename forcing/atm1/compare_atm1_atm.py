# -*- coding: utf-8 -*-
"""
Code to compare the results of atm1 to atm, to make sure
we are doing the variable conversion and interpolation well enough.

"""

gridname='cas6'
tag='v3'
date_string = '2017.04.20'

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

import atm_fun as afun
from importlib import reload
reload(afun)

# get files

grid_fn = Ldir['grid'] + 'grid.nc'
dsg = nc.Dataset(grid_fn)

Lon = dsg['lon_rho'][:]
Lat = dsg['lat_rho'][:]

dsg.close()

indir = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/atm1/'
indiro = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/atm/'
    
#%% plotting
aa = [-130, -122, 42, 52]
lim_dict = dict(zip(afun.outvar_list, afun.lim_list))

plt.close('all')

#outvar_list = ['rain'] # testing
outvar_list = afun.outvar_list

for vn in outvar_list:
    
    ds = nc.Dataset(indir + vn + '.nc')
    dso = nc.Dataset(indiro + vn + '.nc')
    
    vmin = lim_dict[vn][0]
    vmax = lim_dict[vn][1]

    # extract variables to plot
    forecast_hour = 20
    vv = ds[vn][forecast_hour, :, :]
    vvo = dso[vn][forecast_hour, :, :]
    
    ds.close()
    dso.close()
        
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,7), squeeze=False)

    ax = axes[0,0]
    cs = ax.pcolormesh(Lon, Lat, vvo, cmap='rainbow', vmin=vmin, vmax=vmax)
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('Hour = ' + str(forecast_hour) + ', atm')
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
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('atm1')
    ax.axis(aa)

    ax = axes[0,2]
    dvv = vv - vvo
    dv = np.max(np.abs([np.max(dvv), np.min(dvv)]))
    cs = ax.pcolormesh(Lon, Lat, dvv, cmap='bwr', vmin=-dv, vmax=dv)
    fig.colorbar(cs, ax=ax, orientation='horizontal')
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.set_title('atm1 - atm')
    ax.axis(aa)

plt.show()

