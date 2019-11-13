"""
Plot the results of the flux age engine.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zfun

import tef_fun
import flux_fun

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load the DataFrame of results of flux_engine.py
infile = Lfun.choose_item(indir, tag='cc_', itext='Choose flux engine output file:')
cc = pd.read_pickle(indir + infile)

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(indir + 'volumes.p')


plt.close('all')
fig = plt.figure(figsize=(13,8))

cc.loc[:,'age'] = cc.loc[:,'ca']/cc.loc[:,'c']

ax_counter = 1
for ch in flux_fun.channel_dict.keys():
    
    if ax_counter == 1:
        ax = fig.add_subplot(2,1,ax_counter)
    else:
        ax = fig.add_subplot(2,3,ax_counter+2)
        
    seg_list = flux_fun.seg_dict[ch]

    vs = [s + '_s' for s in seg_list]
    vf = [s + '_f' for s in seg_list]
        
    dist = flux_fun.make_dist(v_df.loc[seg_list,'lon'],v_df.loc[seg_list,'lat'])

    # plots of layer tracer concentration or age
    
    if True:
        ax.plot(dist, cc.loc[vs,'c'].values,'-*r')
        ax.plot(dist, cc.loc[vf,'c'].values,'-*r', alpha=.5)
        for ii in range(len(dist)):
            ax.text(dist[ii], cc.loc[vs[ii],'c'], seg_list[ii], color='r')
        ax.set_ylim(bottom=0)
    else:
        # values from the flux_engine
        ax.plot(dist, cc.loc[vs,'age'].values*365,'-*r')
        ax.plot(dist, cc.loc[vf,'age'].values*365,'-*r', alpha=.5)
        for ii in range(len(dist)):
            ax.text(dist[ii], cc.loc[vs[ii],'age']*365, seg_list[ii], color='r')
        ax.set_ylim(0, 700)
        
    if ax_counter == 1:
        ax.set_xlim(-10,410)
    else:
        ax.set_xlim(-10,180)
    
    ax.set_title(ch)
    ax.grid(True)
    ax_counter += 1

plt.show()


