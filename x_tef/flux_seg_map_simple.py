"""
Plot a map of the segments, showing just the key areas.

"""

# imports
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap as lsc
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
from importlib import reload
reload(flux_fun)

# select the indir
indir0 = Ldir['LOo'] + 'tef/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'
# and set the outdir
outdir = indir0 + 'misc_figs_' + Ldir['gridname'] + '/'
Lfun.make_dir(outdir)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# load data
bathy_dict = pickle.load(open(voldir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(voldir + 'ji_dict.p', 'rb'))

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
# if testing:
#     ch_list = [ch_list[0]]

# colors
clist = flux_fun.clist

jj = 0
for ch in ch_list:
    
    cmap = lsc.from_list('pm',[clist[jj],[0,0,0]])

    seg_list = flux_fun.short_seg_dict[ch]
    hh = np.zeros(h.shape)
    for seg_name in seg_list:
        full_ji_list = ji_dict[seg_name]
        for ji in full_ji_list: 
            hh[ji]=1
    H = hh==0
    hm = np.ma.masked_where(H,h)
    vmin=0; vmax=300
    if jj == 0:
        ax.pcolormesh(xp, yp, hm[1:-1,1:-1], cmap=cmap, vmin=vmin, vmax=vmax)
        pfun.dar(ax)
        ax.axis(aa0)
    else:
        ax.pcolormesh(xp, yp, hm[1:-1,1:-1], cmap=cmap, vmin=vmin, vmax=vmax)
    ax.text(.05, .3-jj*.05, ch, color=clist[jj],
        transform=ax.transAxes, fontsize=18, fontweight='bold')
    jj += 1

# add sections and segments
for sn in sect_df.index:
    x0, x1, y0, y1, landward = sect_df.loc[sn,:]
    xx = (x0+x1)/2; yy = (y0+y1)/2
    ax.plot([x0,x1], [y0,y1], '-', color='brown', linewidth=2)
    
ax.set_axis_off()
fig.tight_layout()
plt.show()

if testing == False:
    fig.savefig(outdir + 'seg_map_simple.png')
    




