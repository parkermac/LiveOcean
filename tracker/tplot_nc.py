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
my_npt = int(input('-- Experiment number -- '))
indir = indir_list[my_npt] + '/'
#
rel_list = [rel for rel in os.listdir(indir0 + indir) if 'release' in rel]
rel_list.sort()
Nrl = len(rel_list)
print('\n%s\n' % '** Choose Release file to plot **')
for nrl in range(Nrl):
    print(str(nrl) + ': ' + rel_list[nrl])
my_nrl = int(input('-- Release number -- '))
rel = rel_list[my_nrl]

dsr = nc4.Dataset(indir0 + indir + rel)
dsg = nc4.Dataset(indir0 + indir + 'grid.nc')
    
NT, NP = dsr['lon'].shape

dt_list = [Lfun.modtime_to_datetime(ot) for ot in dsr['ot'][:]]

lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]

aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]

# PLOTTING

plt.close('all')
fig = plt.figure(figsize=(8,8))

ax = fig.add_subplot(111)
pfun.add_coast(ax)

#pfun.add_bathy_contours(ax, G, txt=True)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
# ax.text(.06, .04, ' '.join(p.split('_')),
#     verticalalignment='bottom', transform=ax.transAxes,
#     rotation='vertical')

# add the tracks
ax.plot(dsr['lon'][:], dsr['lat'][:], '-k', alpha = 0.1)
beach_mask = (dsr['u'][-1,:] == 0) & (dsr['v'][-1,:] == 0)
ax.plot(dsr['lon'][:,beach_mask], dsr['lat'][:,beach_mask], '-r', linewidth=1)
ax.plot(dsr['lon'][0,beach_mask],dsr['lat'][0,beach_mask],'or',
        markersize=5, alpha = .4)

plt.show()
#pfun.topfig()
