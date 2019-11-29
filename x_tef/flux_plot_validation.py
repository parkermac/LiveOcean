"""
Plot the results of the flux engine.

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
import flux_fun

# select the indir
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load the DataFrame of TEF transports and salinities, created by flux_make_two_layer.py
df_2 = pd.read_pickle(indir + 'two_layer.p')
# index is ['jdf1', 'jdf2', 'jdf3',...
# columns are ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat']

# load the segment definitions, which is a dictionary defining all the
# sections and rivers surrounding each segment
segs = flux_fun.segs
# this is a dict of dicts, with entries like:
#   {'J1':
#       {'S': [], 'N': [], 'W': ['jdf1'], 'E': ['jdf2'], 'R': ['sanjuan', 'hoko']},...

# load the DataFrame of results of flux_engine.py
cc = pd.read_pickle(indir + 'cc_ocean_salt.p')
# index is ['J1_s', 'J1_f',... = (*)
# columns are ['c', 'v', 'netq']

    
plt.close('all')
fig = plt.figure(figsize=(13,8))

ax_counter = 1
for ch in flux_fun.channel_dict.keys():
    
    if ax_counter == 1:
        ax = fig.add_subplot(2,1,ax_counter)
    else:
        ax = fig.add_subplot(2,3,ax_counter+2)
        
    sect_list = flux_fun.channel_dict[ch]
    seg_list = flux_fun.seg_dict[ch]
    
    # get target salinities (the actual TEF values, at segment boundaries)
    Ns = len(sect_list)
    sa_s = np.zeros(Ns)
    sa_f = np.zeros(Ns)
    sa_x = np.zeros(Ns)
    sa_y = np.zeros(Ns)
    counter = 0
    for sect in sect_list:
        q_s, q_f, f_s, f_f, s_s, s_f, lon, lat = df_2.loc[sect,:]
        sa_s[counter] = s_s
        sa_f[counter] = s_f
        sa_x[counter] = lon
        sa_y[counter] = lat
        counter += 1
    # and associated distance vector
    sa_dist = flux_fun.make_dist(sa_x, sa_y)
    
    vs = [s + '_s' for s in seg_list]
    vf = [s + '_f' for s in seg_list]

    # alternate definition of dist that goes with segments
    if ch in ['Admiralty Inlet to South Sound', 'Hood Canal']:
        dist = sa_dist.copy()
        ddd = np.diff(sa_dist)/2
        dist[:-1] += ddd
        dist[-1] += ddd[-1]
        dist = np.append(0,dist)
    elif ch == 'Whidbey Basin':
        dist = sa_dist[:-1].copy()
        ddd = np.diff(sa_dist)/2
        dist += ddd
        dist = np.append(0,dist)
    else:
        dist = sa_dist[:-1].copy()
        ddd = np.diff(sa_dist)/2
        dist += ddd

    if True:
        # plots of layer salinities
        
        # values from the flux_engine
        ax.plot(dist, cc.loc[vs,'c'].values,'-*r')
        ax.plot(dist, cc.loc[vf,'c'].values,'-*r', alpha=.5)
        for ii in range(len(dist)):
            ax.text(dist[ii], cc.loc[vs[ii],'c'], seg_list[ii], color='r')
        # TEF target values
        ax.plot(sa_dist, sa_s, '-ob')
        ax.plot(sa_dist, sa_f, '-ob', alpha=.5)
        for ii in range(len(sa_dist)):
            ax.text(sa_dist[ii], sa_s[ii], sect_list[ii], color='b')

        if ax_counter == 1:
            ax.set_xlim(-10,410)
        else:
            ax.set_xlim(-10,180)

        if ax_counter == 1:
            ax.set_ylim(28,34)
        elif ax_counter == 2:
            ax.set_ylim(28,32)
        elif ax_counter == 3:
            ax.set_ylim(24,32)
        elif ax_counter == 4:
            ax.set_ylim(22,32)
            

    else:
        # plots of sbar and sprime
        ax2 = ax.twinx()
        
        # values from the flux_engine
        sbot = cc.loc[vs,'c'].values
        stop = cc.loc[vf,'c'].values
        sbar = (sbot+stop)/2
        sprime = sbot - stop
        
        ax.plot(dist, sbar,'-*r')
        ax2.plot(dist, sprime,'-*m')
        for ii in range(len(dist)):
            ax.text(dist[ii], sbar[ii], seg_list[ii], color='r')
        # TEF target values
        sa_sbar = (sa_s+sa_f)/2
        sa_sprime = sa_s - sa_f

        ax.plot(sa_dist, sa_sbar, '-ob')
        ax2.plot(sa_dist, sa_sprime, '-og')
        for ii in range(len(sa_dist)):
            ax.text(sa_dist[ii], sa_sbar[ii], sect_list[ii], color='b')
        
        ax.set_xlim(-10,410)
        ax.set_ylim(22,34)
        
        ax2.set_ylim(0,4)
        
        
    if ax_counter == 4:
        ax.text(.9,.1,'Target',color='b',fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right')
        ax.text(.9,.2,'Model',color='r',fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right')

    ax.text(.05,.05,ch,transform=ax.transAxes)
    
    
    ax.grid(True)
    
    ax_counter += 1

plt.show()


