"""
Code to plot the volumes of the Salish Sea in a graphically compelling way
that we could use for movies of the flux_engine results.

"""

# imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')

import flux_fun

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

outdir = indir0 + 'misc_figs/'

# load a DataFrame of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# vectors of the starting locations of the channels on the plot
x00_list = [0, 17, 8, 19]
y00_list = [0, -9.25, -12, -6]

plt.close('all')
fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111)

ch_list = list(flux_fun.short_seg_dict.keys())

xy = {}

jj = 0
for ch in ch_list:
    seg_list = flux_fun.short_seg_dict[ch].copy()
    
    if ch in ['Hood Canal', 'Whidbey Basin']:
        seg_list.reverse()
            
    # make vectors of volume
    vs = v_df.loc[seg_list,'volume m3'].to_numpy()
    hs = vs**(1/3)
    hs = hs/1e3

    x00 = x00_list[jj]
    y00 = y00_list[jj]
    
    dist = np.cumsum(hs)
    dist = np.append(0,dist)
    ii = 0
    for seg in seg_list:
        x0 = x00 + dist[ii]
        x1 = x00 + dist[ii+1]
        y0 = y00
        y1 = y00 - hs[ii]
        dy = y0-y1
        # bottom layer
        fr = .2
        ax.fill([x0,x1,x1,x0],[y0-fr*dy,y0-fr*dy,y1,y1],color=flux_fun.clist[jj], alpha=.5)
        # top layer
        ax.fill([x0,x1,x1,x0],[y0,y0,y1+(1-fr)*dy,y1+(1-fr)*dy],color=flux_fun.clist[jj], alpha=.5)
        ax.text((x0+x1)/2,y0+.2,seg, horizontalalignment='center', fontsize=10)
        ii += 1
        # save position of center of cell
        xy[seg] = ((x0+x1)/2, (y0+y1)/2)
        
    jj += 1
    
# plot connecting lines
lw = 3
al = .3
ax.plot([xy['J4'][0],xy['A1'][0]], [xy['J4'][1],xy['A1'][1]], '-ok', linewidth=lw, alpha=al)
ax.plot([xy['H1'][0],xy['A3'][0]], [xy['H1'][1],xy['A3'][1]], '-ok', linewidth=lw, alpha=al)
ax.plot([xy['M1'][0],xy['W1'][0]], [xy['M1'][1],xy['W1'][1]], '-ok', linewidth=lw, alpha=al)
ax.plot([xy['W4'][0],xy['J4'][0]], [xy['W4'][1],xy['J4'][1]], '-ok', linewidth=lw-2, alpha=al)
    
ax.set_axis_off()
plt.show()

fig.savefig(outdir + 'volume_plot.png')


