#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 14:22:55 2016

@author: PM5

Plot results of tracker_1, using NetCDF output.
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matplotlib.pyplot as plt

plp = os.path.abspath('../plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import netCDF4 as nc4
import numpy as np

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tracks/'

# choose the file to plot
indir_list_raw = os.listdir(indir0)
indir_list = []
for d in indir_list_raw:
    if os.path.isdir(indir0 + d):
        indir_list.append(d)
indir_list.sort()

Npt = len(indir_list)
#
print('\n%s\n' % '** Choose Experiment to plot **')
for npt in range(Npt):
    print(str(npt) + ': ' + indir_list[npt])
my_npt = input('-- Experiment number (return = 0) --')
if len(my_npt)==0:
    my_npt = 0
indir = indir_list[int(my_npt)] + '/'
#
rel_list = [rel for rel in os.listdir(indir0 + indir) if 'release' in rel]
rel_list.sort()
Nrl = len(rel_list)
print('\n%s\n' % '** Choose Release file to plot **')
for nrl in range(Nrl):
    print(str(nrl) + ': ' + rel_list[nrl])
my_nrl = input('-- Release number (return = 0) -- ')
if len(my_nrl)==0:
    my_nrl = 0
rel = rel_list[int(my_nrl)]

dsr = nc4.Dataset(indir0 + indir + rel)
dsg = nc4.Dataset(indir0 + indir + 'grid.nc')
    
NT, NP = dsr['lon'].shape

dt_list = [Lfun.modtime_to_datetime(ot) for ot in dsr['ot'][:]]

lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]

# PLOTTING

#aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
pad = .02
aa = [dsr['lon'][:].min()-pad, dsr['lon'][:].max()+pad,
    dsr['lat'][:].min()-pad, dsr['lat'][:].max()+pad]

plt.close('all')
fig = plt.figure(figsize=(8,8))

ax = fig.add_subplot(111)
h = dsg['h'][:]
mask = dsg['mask_rho'][:]
zm = -np.ma.masked_where(mask==0, h)
plt.pcolormesh(lonp, latp, zm[1:-1, 1:-1], vmin=-100, vmax=0, cmap='rainbow')
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# add the tracks (packed [time, particle])
ax.plot(dsr['lon'][:], dsr['lat'][:], '-*k')
ax.plot(dsr['lon'][0,:], dsr['lat'][0,:], 'og')
ax.plot(dsr['lon'][-1,:], dsr['lat'][-1,:], 'or')

plt.show()
