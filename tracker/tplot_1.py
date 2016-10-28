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

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
inname = m_dict[my_npt]

import pickle  # python 3
P, G, S, PLdir = pickle.load( open( indir + inname, 'rb' ) )

NT, NP = P['lon'].shape

lonp = G['lon_psi']
latp = G['lat_psi']

if 'greencrab' in inname:
    aa = [-125.5, -122, 47, 49.5]
else:
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
ax.text(.06, .04, ' '.join(inname.split('_')),
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
