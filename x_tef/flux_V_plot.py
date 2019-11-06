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

# load a DataFrame of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# vectors of the starting locations of the channels on the plot
x00_list = [0, 0, 30, 45]
y00_list = [0, -10, -10, -10]

plt.close('all')
fig = plt.figure(figsize=(16, 5))
ax = fig.add_subplot(111)

ch_list = list(flux_fun.seg_dict.keys())

jj = 0
for ch in ch_list:
    seg_list = flux_fun.seg_dict[ch]
    
    # make vectors of volume
    vs = v_df.loc[seg_list,'volume m3'].values
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
        ax.fill([x0,x1,x1,x0],[y0,y0,y1,y1],color=flux_fun.clist[jj], alpha=.5)
        ax.text((x0+x1)/2,y0+.2,seg, horizontalalignment='center')
        ii += 1
    jj += 1
    
ax.set_axis_off()
plt.show()


