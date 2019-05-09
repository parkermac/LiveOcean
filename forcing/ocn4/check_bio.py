# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 07:14:06 2016

@author: PM5

Code to look at the  biogeochemical fields in the ocn files on a given day.

Only set up to work on mac.

"""

gridname='cas6'
tag='v2'
date_string = '2016.12.15'

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

indir = Ldir['LOo'] + Ldir['gtag'] +'/f' + date_string + '/ocn4/'
in_fn = (indir + 'ocean_ini.nc')
ds = nc.Dataset(in_fn)

#%% plotting

vn_dict =  {'salt' : '',
            'temp' : 'deg C',
            'NO3':'MicroMolar',
           'phytoplankton':'MicroMolar N',
           #'zooplankton':'MicroMolar N',
           'detritus':'MicroMolar N',
           #'Ldetritus':'MicroMolar N',
           #'CaCO3':'MicroMolar C',
           'oxygen':'MicroMolar O',
           'alkalinity':'MicroMolar',
           'TIC':'MicroMolar C'}

plt.close()
fig = plt.figure(figsize=(13,7))
count = 1
NR = 2
NC = 4
for vn in vn_dict.keys():
    if vn == 'CaCO3':
        pass
    else:
        ax = fig.add_subplot(NR, NC, count)
        v = ds[vn][0, -1, :, :]
        v[maskr==0] = np.nan
        cs = ax.pcolormesh(lonp, latp, v[1:-1,1:-1], cmap='rainbow')
        ax.set_xticklabels('')
        ax.set_yticklabels('')
        fig.colorbar(cs, ax=ax)
        pfun.dar(ax)
        ax.text(.95, .05, vn,
            horizontalalignment='right', transform=ax.transAxes,
            fontweight='bold')
        if vn == 'phytoplankton':
            ax.text(.05, .9, Ldir['gtag'],
                transform=ax.transAxes, fontweight='bold')
            ax.text(.05, .8, date_string,
                transform=ax.transAxes, fontweight='bold')
        count += 1
plt.show()

