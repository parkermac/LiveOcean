"""
Run the flux engine!
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
import flux_fun

check_convergence = True

# select input/output location
indir0 = '/Users/pm7/Documents/LiveOcean_output/tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# segment definitions, assembled by looking at the figure
# created by plot_thalweg_mean.py
segs = flux_fun.segs

# load DateFrames of transport and volume
q_df = pd.read_pickle(indir + 'q_df.p')
v_df = pd.read_pickle(indir + 'volumes.p')

# create Series of two-layer volumes
V = pd.Series(index=q_df.index)
for seg_name in v_df.index:
    V[seg_name+'_s'] = 0.8 * v_df.loc[seg_name,'volume m3']
    V[seg_name+'_f'] = 0.2 * v_df.loc[seg_name,'volume m3']
    
# specify the forcing array
#
source = 'ocean_salt'
# Valid Choices:
#
# Dye from the ocean meant to reproduce the actual mean salinity
# - ocean_salt
#
# Dye coming in with concentration = 1:
# - ocean
# - river_fraser
# - river_deschutes
# - river_skagit
#
# Initial condition, dye concentration = 1 in a basin
# - ic_hood_canal
#

f = pd.DataFrame(0, index=q_df.index, columns=q_df.columns)

for seg_name in f.index:
    if source == 'ocean_salt':
        if 'J1' in seg_name:
            f.loc[seg_name,'ocean_s'] = 33.1
        elif 'G6' in seg_name:
            f.loc[seg_name,'ocean_s'] = 30.5
    elif source == 'ocean':
        if 'J1' in seg_name:
            f.loc[seg_name,'ocean_s'] = 1
        elif 'G6' in seg_name:
            f.loc[seg_name,'ocean_s'] = 1
    elif source == 'river_fraser':
        if 'G3' in seg_name:
            f.loc[seg_name,'river_f'] = 1
    elif source == 'river_deschutes':
        if 'S4' in seg_name:
            f.loc[seg_name,'river_f'] = 1
    elif source == 'river_skagit':
        if 'W4' in seg_name:
            f.loc[seg_name,'river_f'] = 1

# we will do the integration just with numpy arrays
NR = len(q_df.index)
NC = len(q_df.columns)
vv = V.values
ivv = 1/vv

c = np.zeros(NR)
ca = np.zeros(NR) 
q = q_df.values
ff = np.zeros((NR,NC))
ffa = np.zeros((NR,NC))

if 'ocean' in source:
    ff[:,0] = f.loc[:,'ocean_s'].values # column 0 is the ocean inflow
elif 'river' in source:
    ff[:,3] = f.loc[:,'river_f'].values # column 3 is the river inflow

dt = 3e3 # time step (seconds)
nyears = 6
NT = int(nyears*365*86400/dt) # number of time steps

if check_convergence:
    c_check = np.nan + np.ones((int(NT/100)+1, NR))
    t_check = np.nan + np.ones(int(NT/100)+1)
    
sinking = False
dz = 1e-4 # a parameter to control sinking rate

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
        
    qffa = (q*ffa).sum(axis=1)
    ca = ca + dt*ivv*qffa + dt*c/(365*86400)
    ffa[:,4:] = np.tile(ca,(NR,1))
    
    if check_convergence and np.mod(ii,100)==0 :
        c_check[int(ii/100), :] = c.copy()
        t_check[int(ii/100)] = dt * ii
    
if check_convergence:
    plt.close('all')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.plot(t_check/(365*86400), c_check, '-')
    ax.set_xlabel('Time (years)')
    ax.set_xlim(0,nyears)
    plt.show()

# create and fill the output DataFrame
cc = pd.DataFrame(index=q_df.index,columns=['c', 'ca'])
cc['c'] = c
cc['ca'] = ca
        
cc.to_pickle(indir + 'cc_' + source + '.p')
    
    