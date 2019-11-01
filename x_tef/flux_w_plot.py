"""
Code to plot the spatial distribution of the vertical velocities that
result from flux_get_A.py.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zfun

import tef_fun
from importlib import reload
reload(tef_fun)

import flux_fun
reload(flux_fun)

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# this is the big DataFrame created by flux_get_A.py
q_df = pd.read_pickle(indir + 'q_df.p')
# index is ['J1_s', 'J1_f', 'J2_s',... = (*)
# columns are ['ocean_s', 'ocean_f', 'river_s', 'river_f', 'J1_s', 'J1_f', 'J2_s',...

# make lists of the various segment sequences
ssJ = ['J'+str(s) for s in range(1,5)]
ssM = ['M'+str(s) for s in range(1,10)]
ssMH = ['M'+str(s) for s in range(1,4)] # to Hood Canal
ssMW = ['M'+str(s) for s in range(1,5)] # to Whidbey Basin
ssS = ['S'+str(s) for s in range(1,7)]
ssG = ['G'+str(s) for s in range(1,7)]
ssW = ['W'+str(s) for s in range(1,5)]
ssH = ['H'+str(s) for s in range(1,9)]

def make_dist(x,y):
    NS = len(x)
    xs = np.zeros(NS)
    ys = np.zeros(NS)
    xs, ys = zfun.ll2xy(x, y, x[0], y[0])
    dx = np.diff(xs)
    dy = np.diff(ys)
    dd = (dx**2 + dy**2)**.5 # not clear why np.sqrt throws an error
    dist = np.zeros(NS)
    dist[1:] = np.cumsum(dd/1000)
    return dist
    
seg_dict = {'JdF to South Sound': ssJ + ssM + ssS,
            'JdF to Strait of Georgia': ssJ + ssG,
            'JdF to Hood Canal': ssJ + ssMH + ssH,
            'JdF to Whidbey Basin': ssJ + ssMW + ssW}

#plt.close('all')
fig = plt.figure(figsize=(13,8))
lw = 3
fs = 16
abc = 'abcd'

ax_counter = 1
for ch in seg_dict.keys():
    
    ax = fig.add_subplot(2,2,ax_counter)
    ax.tick_params(labelsize=fs-2) # tick labels

    seg_list = seg_dict[ch]
    dist = make_dist(v.loc[seg_list,'lon'],v.loc[seg_list,'lat'])
        
    # make vectors of vertical velocity on the segments of this channel
    w_up = np.nan * np.ones(len(seg_list))
    w_down = np.nan * np.ones(len(seg_list))
    seg_counter = 0
    for seg in seg_list:
        w_up[seg_counter] = q_df.loc[seg+'_f',seg+'_s'] / v.loc[seg,'area m2']
        w_down[seg_counter] = q_df.loc[seg+'_s',seg+'_f'] / v.loc[seg,'area m2']
        seg_counter += 1 
    
    # plotting
    ax.plot(dist, w_up*1e3,'-*r', linewidth=lw)
    ax.plot(dist, w_down*1e3,'-*b', linewidth=lw)
    # for ii in range(len(dist)):
    #     ax.text(dist[ii], w_up[ii]*1e3, seg_list[ii], color='r')

    # formatting and labels
    ax.set_xlim(-10,500)
    if ax_counter == 1:
        ax.set_ylim(0, .7)
        ax.text(.05,.75,'Note different y-scale',transform=ax.transAxes, fontsize=fs-2, style='italic')
    else:
        ax.set_ylim(0,.25)
    if ax_counter in [1,3]:
        ax.set_ylabel('Vertical Velocity (mm/s)', fontsize=fs)
    if ax_counter in [3,4]:
        ax.set_xlabel('Distance from JdF (km)', fontsize=fs)
    if ax_counter in [2,4]:
        pass
        #ax.set_yticklabels([])
    if ax_counter in [1,2]:
        ax.set_xticklabels([])
    if ax_counter == 4:
        ax.text(.9,.2,'Upwards',color='r',fontweight='bold',transform=ax.transAxes, fontsize=fs, horizontalalignment='right')
        ax.text(.9,.1,'Downwards',color='b',fontweight='bold',transform=ax.transAxes, fontsize=fs, horizontalalignment='right')
    ax.text(.05, .9, '('+abc[ax_counter-1]+') '+ch , fontsize=fs, transform=ax.transAxes)
    ax.grid(True)
    ax_counter += 1

plt.show()


