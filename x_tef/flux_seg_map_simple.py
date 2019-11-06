"""
Plot a map of the segments, showing just the key areas.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun

gridname = 'cas6'
tag = 'v3'
Ldir = Lfun.Lstart(gridname, tag)

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load data
bathy_dict = pickle.load(open(indir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))

h = bathy_dict['h']
xp = bathy_dict['xp']
yp = bathy_dict['yp']

testing = False

# plotting
plt.close('all')
fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(111)
aa0 = [-125.5, -122, 46.7, 50.4]

ch_list = list(flux_fun.seg_dict.keys())
if testing:
    ch_list = [ch_list[0]]

# colors
clist = flux_fun.clist
cmap_list = []
for c in clist:
    cmap_list.append(c.title() + 's')

jj = 0
for ch in ch_list:
    seg_list = flux_fun.seg_dict[ch]
    hh = np.zeros(h.shape)
    for seg_name in seg_list:
        full_ji_list = ji_dict[seg_name]
        for ji in full_ji_list: 
            hh[ji]=1
    H = hh==0
    hm = np.ma.masked_where(H,h)
    vmin=-100; vmax=200
    if jj == 0:
        ax.pcolormesh(xp, yp, hm[1:-1,1:-1], cmap=cmap_list[jj], vmin=vmin, vmax=vmax)
        pfun.dar(ax)
        ax.axis(aa0)
    else:
        ax.pcolormesh(xp, yp, hm[1:-1,1:-1], cmap=cmap_list[jj], vmin=vmin, vmax=vmax)
    ax.text(.05, .3-jj*.05, ch, color=clist[jj],
        transform=ax.transAxes, fontsize=16)
    jj += 1
    
ax.set_axis_off()
fig.tight_layout()
plt.show()
    




