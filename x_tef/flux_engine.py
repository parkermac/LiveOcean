"""
Find the volume of each flux segment
"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import netCDF4 as nc

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zrfun

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
from importlib import reload
reload(tef_fun)

import flux_fun
reload(flux_fun)

verbose = False

# select output location
if False:
    indir00 = Ldir['LOo'] + 'tef/'
    # choose the tef extraction to work with
    item = Lfun.choose_item(indir00)
    indir = indir0 + item + '/'
else:
    indir0 = '/Users/pm7/Documents/LiveOcean_output/tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

df = pd.read_pickle(indir + 'two_layer.p')
#q_s, q_f, f_s, f_f, s_s, s_f, lon, lat = df.loc[sect,:]
# segment definitions, assembled by looking at the figure
# created by plot_thalweg_mean.py
segs = flux_fun.segs


q_df = pd.read_pickle(indir + 'q_df.p')
volumes = pd.read_pickle(indir + 'volumes.p')

V = pd.Series(index=q_df.index)
for seg_name in volumes.index:
    V[seg_name+'_s'] = 0.8 * volumes.loc[seg_name,'volume m3']
    V[seg_name+'_f'] = 0.2 * volumes.loc[seg_name,'volume m3']
    
# the forcing array
f = pd.DataFrame(0, index=q_df.index, columns=q_df.columns)
for seg_name in f.index:
    if 'J1' in seg_name:
        f.loc[seg_name,'ocean_s'] = 33.1
    elif 'G6' in seg_name:
        f.loc[seg_name,'ocean_s'] = 30.5

NR = len(q_df.index)
NC = len(q_df.columns)

# we will do the integration just with numpy arrays
vv = V.values
ivv = 1/vv
c = np.zeros(NR)
q = q_df.values
ff = np.zeros((NR,NC))
ff[:,0] = f.loc[:,'ocean_s'].values

dt = 3e3
for ii in range(30000):
    #print(ii)
    
    qff = (q*ff).sum(axis=1)
    c = c + dt*ivv*qff
    ff[:,4:] = np.tile(c,(NR,1))
    

cc = pd.DataFrame(index=q_df.index,columns=['c', 'v', 'netq'])
cc['c'] = c
cc['v'] = vv
cc['netq'] = q.sum(axis=1)

if verbose:
    for seg_name in cc.index:
        print('%5s: %10.2f %10.2f %10.2f' % (seg_name,
            cc.loc[seg_name,'c'], cc.loc[seg_name,'v']/1e9, cc.loc[seg_name,'netq']/1e3))
        
vs = (['J'+str(s)+'_s' for s in range(1,5)] +
        ['M'+str(s)+'_s' for s in range(1,10)]+
        ['S'+str(s)+'_s' for s in range(1,4)])
vf = (['J'+str(s)+'_f' for s in range(1,5)] +
        ['M'+str(s)+'_f' for s in range(1,10)]+
        ['S'+str(s)+'_f' for s in range(1,4)])

# plotting
plt.close('all')
fig = plt.figure(figsize=(14,8))

ax = fig.add_subplot(111)
cc.loc[vs,'c'].plot(ax=ax, color='r')
cc.loc[vf,'c'].plot(ax=ax, color='b')

plt.show()    
    
    