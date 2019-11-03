"""
Run the flux engine to determine the age of different water sources.
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

check_convergence = True

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
source = 'ocean'
f = pd.DataFrame(0, index=q_df.index, columns=q_df.columns)
for seg_name in f.index:
    if source == 'ocean':
        if 'J1' in seg_name:
            f.loc[seg_name,'ocean_s'] = 1
        elif 'G6' in seg_name:
            f.loc[seg_name,'ocean_s'] = 1
    elif source == 'river':
        if 'G3' in seg_name: # G3 = Fraser, S4 = Deschutes, W4 = Skagit
            f.loc[seg_name,'river_f'] = 1

NR = len(q_df.index)
NC = len(q_df.columns)

# we will do the integration just with numpy arrays
vv = V.values
ivv = 1/vv
c = np.zeros(NR)
ca = np.zeros(NR) 
q = q_df.values
ff = np.zeros((NR,NC))
ffa = np.zeros((NR,NC))
if source == 'ocean':
    ff[:,0] = f.loc[:,'ocean_s'].values # column 0 is the ocean inflow
elif source == 'river':
    ff[:,3] = f.loc[:,'river_f'].values # column 3 is the river inflow

dt = 3e3 # time step (seconds)
NT = 2*60000 # number of time steps

if check_convergence:
    c_check = np.nan + np.ones((int(NT/100), NR))
    t_check = np.nan + np.ones(int(NT/100))
    
sinking = True
dz = 1e-4 # a parameter to control sinking rate

# calculate distribution of tracer
for ii in range(NT):
    qff = (q*ff).sum(axis=1)
    c = c + dt*ivv*qff
    
    if sinking == True:
        NC2 = int(len(c)/2)
        for jj in range(NC2):
            c_f = c[2*jj + 1]
            c[2*jj + 1] -= c_f*dz
            c[2*jj] += c_f*dz
            
        
    ff[:,4:] = np.tile(c,(NR,1))
    if check_convergence and np.mod(ii,100)==0 :
        c_check[int(ii/100), :] = c.copy()
        t_check[int(ii/100)] = dt * ii
        
    qffa = (q*ffa).sum(axis=1)
    ca = ca + dt*ivv*qffa + dt*c/(365*86400)
    ffa[:,4:] = np.tile(ca,(NR,1))
    
        
if check_convergence:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t_check/86400, c_check, '-')
    plt.show()
    
cc = pd.DataFrame(index=q_df.index,columns=['c', 'ca', 'v', 'netq'])
cc['c'] = c
cc['ca'] = ca
cc['v'] = vv
cc['netq'] = q.sum(axis=1)

        
cc.to_pickle(indir + 'cc_age.p')
    
    