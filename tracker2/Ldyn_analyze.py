"""
This code loads in the results of Ldyn_gather and analyzes the
along-track momentum balance of the "winners."

Performance:

For testing run as:

"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun

Ldir = Lfun.Lstart()

import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
from time import time
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import Ldyn_functions as Ldf
from importlib import reload
reload(Ldf)

# command line inputs
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-exp_name', type=str, default='EJdF3d_ndiv12_3d_up1')
parser.add_argument('-testing', type=zfun.boolean_string, default=False)
parser.add_argument('-verbose', type=zfun.boolean_string, default=False)
args = parser.parse_args()
exp_name = args.exp_name
testing = args.testing
verbose = args.verbose

t_dir = Ldir['LOo'] + 'tracks2/' + exp_name + '/'
EI = Lfun.csv_to_dict(t_dir + 'exp_info.csv')
t_fn = t_dir + 'release_' + EI['ds_first_day'] + '.nc'

# get the tracker2 output for this release
t_ds = nc.Dataset(t_fn)
NT, NP = t_ds['lon'].shape
# get a list of datetimes
ot_vec = t_ds['ot'][:]
t_dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]
# Gather particle data; packed [time, particle #]
lon = t_ds['lon'][:]
lat = t_ds['lat'][:]
cs = t_ds['cs'][:]
t_ds.close()

# load results of Ldyn_gather.py
u_dict = pickle.load(open(t_dir + 'Ldyn_u_dict.p', 'rb'))
v_dict = pickle.load(open(t_dir + 'Ldyn_v_dict.p', 'rb'))
imask = pickle.load(open(t_dir + 'Ldyn_imask.p', 'rb'))

a_list = ['accel','hadv', 'vadv','cor','prsgrd','prsgrd0','vvisc']

# make sure everything is numeric
for vn in a_list:
    for cn in u_dict['u_'+vn].columns:
        u_dict['u_'+vn][cn] = pd.to_numeric(u_dict['u_'+vn][cn])
    for cn in v_dict['v_'+vn].columns:
        v_dict['v_'+vn][cn] = pd.to_numeric(v_dict['v_'+vn][cn])

# find dx and dy for the particles
wlon = zfun.fillit(lon[:,imask])
wlat = zfun.fillit(lat[:,imask])
wcs = zfun.fillit(cs[:,imask])
x, y = zfun.ll2xy(wlon,wlat, wlon.mean(), wlat.mean())
dx = np.diff(x, axis=0)
dy = np.diff(y, axis=0)
dd = np.sqrt(dx**2 + dy**2)

# rotate into along-track momeuntum budget
u_cor = u_dict['u_cor'].to_numpy()
v_cor = v_dict['v_cor'].to_numpy()
theta = np.arctan2(u_cor, -v_cor) # radians
sin_th = np.sin(theta)
cos_th = np.cos(theta)

s_dict = {}
n_dict = {}
# "s" stand's for along-track, "n" stands for cross-track
for vn in a_list:
    s_dict['s_'+vn] = cos_th*u_dict['u_'+vn] + sin_th*v_dict['v_'+vn]
    n_dict['n_'+vn] = cos_th*v_dict['v_'+vn] - sin_th*u_dict['u_'+vn]
    
# form along-track spatial integrals
si_dict = {}
ni_dict = {}
for vn in a_list:
    si_dict['s_'+vn] = (dd*s_dict['s_'+vn]).cumsum(axis=0)
    ni_dict['n_'+vn] = (dd*n_dict['n_'+vn]).cumsum(axis=0)
    
# form more informative groupings of terms
SI_dict = {}
SI_dict['accel'] = si_dict['s_accel'] - si_dict['s_hadv'] - si_dict['s_vadv'] # LHS
SI_dict['pg'] = si_dict['s_prsgrd'] # RHS
SI_dict['fric'] = si_dict['s_vvisc'] # RHS
SI_dict['sum'] = SI_dict['accel'] - SI_dict['pg'] - SI_dict['fric']

# PLOTTING
plt.close('all')
fs = 16
plt.rc('font', size=fs)

for ii in range(5):
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111)
    for vn in SI_dict.keys():
        SI_dict[vn][imask[ii]].plot(ax=ax, legend=True, label=vn)
    
plt.show()
plt.rcdefaults()

