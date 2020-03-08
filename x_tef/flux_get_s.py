"""
A tool to extract hourly time series of salinity and volume in the segments.

Performance: takes 1.5 sec per save (3.6 hours per year) on my mac,
and about 50 percent faster on perigee.  This relies
on creating i_dict and j_dict of indices used for fancy indexing in the segment loop.
The alternate version takes about 15 times longer.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import netCDF4 as nc
import argparse
from datetime import datetime, timedelta

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun

import tef_fun
import flux_fun

from time import time

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.07.04')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.07.04')
args = parser.parse_args()

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
# get time limits
ds0 = args.date_string0; ds1 = args.date_string1
Ldir['date_string0'] = ds0; Ldir['date_string1'] = ds1
dt0 = datetime.strptime(ds0, '%Y.%m.%d'); dt1 = datetime.strptime(ds1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1
print('Working on:')
print(Ldir['gtagex'] + '_' + ds0 + '_' + ds1 +'\n')

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
NT = len(fn_list)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = G['h']
DA = G['DX'] * G['DY']
DA3 = DA.reshape((1,G['M'],G['L']))

# select input/output location
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DataFrame of volume and associated dicts
v_df = pd.read_pickle(indir + 'volumes.p')
bathy_dict = pickle.load(open(indir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))

seg_list = list(v_df.index)
verbose = False

testing = False
if testing:
    verbose = True
    seg_list = seg_list[-2:]

j_dict = {}; i_dict = {}
for seg_name in seg_list:
    jj = []; ii = []
    ji_list_full = ji_dict[seg_name]
    for ji in ji_list_full:
        jj.append(ji[0])
        ii.append(ji[1])
    jjj = np.array(jj)
    iii = np.array(ii)
    j_dict[seg_name] = jjj
    i_dict[seg_name] = iii

s_df = pd.DataFrame(columns=seg_list)
v_df = pd.DataFrame(columns=seg_list)

for fn in fn_list:
    
    tt0 = time()
            
    print(fn)
        
    ds = nc.Dataset(fn)
    salt = ds['salt'][0,:,:,:]
    zeta = ds['zeta'][0,:,:]
    ot = ds['ocean_time'][:]
    ds.close()
    
    if testing:
        z_w_alt = zrfun.get_z(h, zeta, S, only_w=True)
        dz_alt = np.diff(z_w_alt, axis=0)
        DV_alt = dz_alt * DA3
    
    dt = Lfun.modtime_to_datetime(ot.data[0])
    
    # find the volume and volume-mean salinity
    for seg_name in seg_list:
        
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        
        z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S, only_w=True)
        dz = np.diff(z_w, axis=0)
        DV = dz * DA3[0,jjj,iii]
        volume = DV.sum()
        net_salt = (salt[:,jjj,iii] * DV).sum()
        mean_salt = net_salt/volume
        
        if testing:
            net_salt_alt = 0
            volume_alt = 0
            full_ji_list = ji_dict[seg_name]
            for ji in full_ji_list:
                net_salt_alt += (salt[:,ji[0],ji[1]] * DV_alt[:,ji[0],ji[1]]).sum()
                volume_alt += DV_alt[:,ji[0],ji[1]].sum()
            mean_salt_alt = net_salt_alt / volume_alt
        
        # store results
        s_df.loc[dt, seg_name] = mean_salt
        v_df.loc[dt, seg_name] = volume
        
        if verbose:
            print('%3s: Mean Salinity = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_salt, volume/1e9))
            print('%3s: Mean Sal alt  = %0.4f, Vol alt = %0.4f km3' %
                (seg_name, mean_salt_alt, volume_alt/1e9))
                
    print('  ** took %0.1f sec' % (time()-tt0))

s_out_fn = indir + 'hourly_segment_salinity.p'
v_out_fn = indir + 'hourly_segment_volume.p'
s_df.to_pickle(s_out_fn)
v_df.to_pickle(v_out_fn)
    
        
    
