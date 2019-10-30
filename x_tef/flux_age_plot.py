"""
Plot the results of the flux age engine.

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

# load a Series of the volumes of each segment, created by flux_get_vol.py
v = pd.read_pickle(indir + 'volumes.p')
# index is ['J1', 'J2', 'J3',...
# columns are ['volume m3', 'area m2', 'lon', 'lat']

# this is the big DataFrame created by flux_get_A.py
q_df = pd.read_pickle(indir + 'q_df.p')
# index is ['J1_s', 'J1_f', 'J2_s',... = (*)
# columns are ['ocean_s', 'ocean_f', 'river_s', 'river_f', 'J1_s', 'J1_f', 'J2_s',...

# Make a version of "v" which has entries for the _s (salty) and _f (fresh)
# parts of each segment, to be compatible with "cc" below.  The new version
V = pd.DataFrame(index=q_df.index, columns=['v_s', 'v_f', 'lon', 'lat'])
for seg_name in v.index:
    V.loc[seg_name+'_s','lon'] = v.loc[seg_name,'lon']
    V.loc[seg_name+'_s','lat'] = v.loc[seg_name,'lat']
    V.loc[seg_name+'_f','lon'] = v.loc[seg_name,'lon']
    V.loc[seg_name+'_f','lat'] = v.loc[seg_name,'lat']

# load the DataFrame of results of flux_engine.py
cc = pd.read_pickle(indir + 'cc_age.p')
# index is ['J1_s', 'J1_f',... = (*)
# columns are ['c', 'ca', 'v', 'netq']

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
    

# also cue up a line for the input salinities from the TEF sections
channel_dict = {'JdF to South Sound':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3','ai4',
                'mb1','mb2','mb3','mb4','mb5',
                'tn1','tn2','tn3',
                'ss1','ss2','ss3'],
            'JdF to Strait of Georgia':['jdf1','jdf2','jdf3','jdf4',
                'sji1', 'sji2', 'sog1','sog2','sog3','sog4','sog5'],
            'JdF to Hood Canal':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3',
                'hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
            'JdF to Whidbey Basin':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3', 'ai4',
                'wb1','wb2','wb3','wb4','dp']}
                
seg_dict = {'JdF to South Sound': ssJ + ssM + ssS,
            'JdF to Strait of Georgia': ssJ + ssG,
            'JdF to Hood Canal': ssJ + ssMH + ssH,
            'JdF to Whidbey Basin': ssJ + ssMW + ssW}

#plt.close('all')
fig = plt.figure(figsize=(13,8))

cc.loc[:,'age'] = cc.loc[:,'ca']/cc.loc[:,'c']

ax_counter = 1
for ch in channel_dict.keys():
    
    ax = fig.add_subplot(2,2,ax_counter)    

    seg_list = seg_dict[ch]
    vs = [s + '_s' for s in seg_list]
    vf = [s + '_f' for s in seg_list]
        
    dist = make_dist(v.loc[seg_list,'lon'],v.loc[seg_list,'lat'])

    # plots of layer tracer concentration or age
    
    if False:
        # values from the flux_engine
        ax.plot(dist, cc.loc[vs,'c'].values,'-*r')
        ax.plot(dist, cc.loc[vf,'c'].values,'-*r', alpha=.5)
        for ii in range(len(dist)):
            ax.text(dist[ii], cc.loc[vs[ii],'c'], seg_list[ii], color='r')
        ax.set_ylim(0, 1)
    else:
        # values from the flux_engine
        ax.plot(dist, cc.loc[vs,'age'].values,'-*r')
        ax.plot(dist, cc.loc[vf,'age'].values,'-*r', alpha=.5)
        for ii in range(len(dist)):
            ax.text(dist[ii], cc.loc[vs[ii],'age']*365, seg_list[ii], color='r')
        #ax.set_ylim(0, 365)
    ax.set_xlim(-10,500)
    ax.set_title(ch)
    ax.grid(True)
    ax_counter += 1

plt.show()


