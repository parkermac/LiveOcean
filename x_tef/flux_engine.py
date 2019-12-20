"""
Run the flux engine!
"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
#import netCDF4 as nc
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import zrfun

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import tef_fun
import flux_fun

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-src', '--source', nargs='?', type=str, default='ocean')
parser.add_argument('-sink', '--sinking', default=False, type=boolean_string)
parser.add_argument('-conv', '--check_convergence', default=False, type=boolean_string)
args = parser.parse_args()

source = args.source
sinking = args.sinking
check_convergence = args.check_convergence

print('Running integration with source = ' + source)
# Valid choices for source:
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

# select input/output location
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DateFrames of transport and volume
q_df = pd.read_pickle(indir + 'q_df.p')
v_df = pd.read_pickle(indir + 'volumes.p')

# create Series of two-layer volumes
V = pd.Series(index=q_df.index)
for seg_name in v_df.index:
    V[seg_name+'_s'] = 0.8 * v_df.loc[seg_name,'volume m3']
    V[seg_name+'_f'] = 0.2 * v_df.loc[seg_name,'volume m3']
    
# "f" is a DataFrame organized like q_df but whose entries
# are the forced values of the tracer
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
            f.loc[seg_name,'ocean_s'] = 30.5/33.1
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
    
if 'ic' in source:
    for seg_name in f.index:
        if source == 'ic_hood_canal':
            this_seg_list = ['H'+ str(item) for item in range(3,9)]
            if seg_name[:2] in this_seg_list:
                jj = int(np.argwhere(f.index==seg_name))
                ff[jj,jj+4] = 1
                c[jj] = 1

dt = 3600 #3e3 # time step (seconds)
nyears = 1#6
NT = int(nyears*365*86400/dt) # number of time steps
savedays = 2#10
Nsave = int(savedays*86400/dt)
# Nsave = number of time steps between saves in order to save every "savedays"

# arrays to store time-varying information
c_arr = np.nan + np.ones((int(NT/Nsave)+1, NR))
t_arr = np.nan + np.ones(int(NT/Nsave)+1)
    
dz = .5e-4 #1e-4 # a parameter to control sinking rate

for ii in range(NT):
    
    if np.mod(ii,Nsave)==0 :
        c_arr[int(ii/Nsave), :] = c.copy()
        t_arr[int(ii/Nsave)] = dt * ii
    
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
        
if check_convergence:
    plt.close('all')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.plot(t_arr/(365*86400), c_arr, '-')
    ax.set_xlabel('Time (years)')
    ax.set_xlim(0,nyears)
    plt.show()

# create and fill the output DataFrame
cc = pd.DataFrame(index=q_df.index,columns=['c', 'ca'])
cc['c'] = c
cc['ca'] = ca

# and another one for the time-dependent fields
aa = pd.DataFrame(c_arr, index=t_arr/86400, columns=q_df.index)
# note that the index is time (days)

sink_tag = ''
if sinking:
    sink_tag = '_sinking'
        
cc.to_pickle(indir + 'cc_' + source + sink_tag + '.p')
aa.to_pickle(indir + 'aa_' + source + sink_tag + '.p')
    
    