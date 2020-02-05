"""
Code to make a nice plot for Figure 1 of the first LiveOcean paper.

"""

# imports

#%% setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun, zrfun
import pfun

from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
import matplotlib.pyplot as plt


fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f2019.07.04/ocean_his_0020.nc'
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G = True)
lon = G['lon_psi'][:]
lat = G['lat_psi'][:]

plt.close('all')

lwt = .5
fs = 14

#vn = 'salt'; vmin=20; vmax=33
vn = 'temp'; vmin=8; vmax=20

v = ds[vn][0,-1,1:-1,1:-1]
cmap = 'RdYlBu_r'#'cubehelix'#'ocean'

fig = plt.figure(figsize=(18,8))

ax = fig.add_subplot(141)
ax.pcolormesh(lon,lat,v, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.add_bathy_contours(ax, ds, txt=True)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(ds))
pfun.dar(ax)

ax.set_xticks([-130, -126, -122])
ax.set_yticks([42, 44, 46, 48, 50, 52])
ax.tick_params(labelsize=fs) # tick labels

# add ticks for grid spacing
x = lon[0,::10]
y = lat[::10,0]
for xx in x:
    ax.plot([xx,xx],[42,42.1],'-k', lw=lwt)
for yy in y:
    ax.plot([-122.2, -122],[yy,yy],'-k', lw=lwt)

ax = fig.add_subplot(142)
ax.pcolormesh(lon,lat,v, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.add_coast(ax)
ax.axis([-123.6, -122, 47, 49])
pfun.dar(ax)

ax.set_xticks([-123, -122])
ax.set_yticks([47, 48, 49])
ax.tick_params(labelsize=fs) # tick labels

# add ticks for grid spacing
x = lon[0,::4]
y = lat[::4,0]
for xx in x:
    ax.plot([xx,xx],[47,47.02],'-k', lw=lwt)
for yy in y:
    hh = ax.plot([-122.04, -122],[yy,yy],'-k', lw=lwt)

fig.tight_layout()
        
plt.show()

