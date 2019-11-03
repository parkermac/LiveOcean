"""
Code to plot the volumes of the Salish Sea in a graphically compelling way
that we could use for movies of the flux_engine results.

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
    
            
x00_list = [0, 0, 30, 45]
y00_list = [0, -10, -10, -10]

plt.close('all')
fig = plt.figure(figsize=(16, 5))
lw = 3
fs = 16
abc = 'abcd'

ax = fig.add_subplot(111)

ax_counter = 1

testing = True

ch_list = list(flux_fun.seg_dict.keys())
# if testing:
#     ch_list = [ch_list[0]]

cvec = 'rbgc'
jj = 0
for ch in ch_list:
    
    seg_list = flux_fun.seg_dict[ch]
    # dist = make_dist(v.loc[seg_list,'lon'],v.loc[seg_list,'lat'])
    # # make an array distance for segment edges
    # dd = np.diff(dist)
    # de = np.append(dist, dist[-1] + dd[-1])
    # dde = np.diff(de)
        
    # make vectors of volume
    vs = v.loc[seg_list,'volume m3'].values
    hs = vs**(1/3)
    hs = hs/1e3
    
    x00 = x00_list[jj]
    y00 = y00_list[jj]
    
    # plotting
    ii = 0
    dist = np.cumsum(hs)
    dist = np.append(0,dist)
    
    for seg in seg_list:
        x0 = x00 + dist[ii]
        x1 = x00 + dist[ii+1]
        y0 = y00
        y1 = y00 - hs[ii]
        # # print(y1)
        
        # #ax.plot(dist, vs, '-o', linewidth=3)
        ax.fill([x0,x1,x1,x0],[y0,y0,y1,y1],color=cvec[jj], alpha=.5)

        ax.text((x0+x1)/2,y0+.2,seg, horizontalalignment='center')
        
        ii += 1


    jj += 1
    
ax.set_axis_off()

plt.show()


