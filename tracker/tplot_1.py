#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 14:22:55 2016

@author: PM5

Plot results of tracker.

Need to update with new plotting functions.
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

import pickle
import numpy as np

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tracks/'

# choose the type of plot to make
print('\n%s\n' % '** Choose directory to plot **')
indir_list_raw = os.listdir(indir0)
indir_list = []
for d in indir_list_raw:
    if os.path.isdir(indir0 + d):
        indir_list.append(d)
Npt = len(indir_list)
indir_dict = dict(zip(range(Npt), indir_list))
for npt in range(Npt):
    print(str(npt) + ': ' + indir_list[npt])
my_npt = int(input('-- Input number -- '))
indir = indir_dict[my_npt]

p_list = os.listdir(indir0 + indir)

p_list.sort()

counter = 0
P = dict()
for p in p_list:

    if counter == 0:
        Pp, G, S, PLdir = pickle.load( open( indir0 + indir + '/' + p, 'rb' ) )
        for k in Pp.keys():
            P[k] = Pp[k]
    else:
        Pp, PLdir = pickle.load( open( indir0 + indir + '/' + p, 'rb' ) )
        for k in Pp.keys():
            if k == 'ot':
                P[k] = np.concatenate((P[k], Pp[k][1:]), axis=0)
            else:
                P[k] = np.concatenate((P[k], Pp[k][1:,:]), axis=0)
    counter += 1
    
NT, NP = P['lon'].shape

dt_list = [Lfun.modtime_to_datetime(ot) for ot in P['ot']]

lonp = G['lon_psi']
latp = G['lat_psi']

aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]

# PLOTTING

plt.close()
fig = plt.figure(figsize=(12,12))

ax = fig.add_subplot(111)
pfun.add_coast(ax)

pfun.add_bathy_contours(ax, G, txt=True)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.text(.06, .04, ' '.join(p.split('_')),
    verticalalignment='bottom', transform=ax.transAxes,
    rotation='vertical')

# add the tracks
ax.plot(P['lon'], P['lat'], '-k', alpha = 0.1)
beach_mask = P['u'][-1,:] == 0
ax.plot(P['lon'][:,beach_mask], P['lat'][:,beach_mask], '-r', linewidth=1)
ax.plot(P['lon'][0,beach_mask],P['lat'][0,beach_mask],'or',
        markersize=5, alpha = .4)

plt.show()
pfun.topfig()
