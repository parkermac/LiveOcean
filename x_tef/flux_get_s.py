"""
A tool to extract time series of salinity statistics in the segment volumes.

We do this with zero SSH and just once per day, because we are only concerned
with low-frequency salinity variation.

A variant of this code could be used to get hourly values with SSH variation.
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

z_w = zrfun.get_z(h, 0*h, S, only_w=True)
dz = np.diff(z_w, axis=0)
DV = dz * DA3

# select input/output location
indir0 = '/Users/pm7/Documents/LiveOcean_output/tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DateFrames of transport and volume
v_df = pd.read_pickle(indir + 'volumes.p')

# load data
bathy_dict = pickle.load(open(indir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))

seg_list = list(v_df.index)

dti = pd.date_range(dt0, dt1, freq='D')

s_df = pd.DataFrame(index=dti, columns=seg_list)

verbose = False

for dt in s_df.index:
    
    tt0 = time()
        
    ymd = datetime.strftime(dt, '%Y.%m.%d')
    #hrs = ('0000' + str(dt.hour + 1))[-4:]
    hrs = '0013'
    fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
        '/f' + ymd + '/ocean_his_' + hrs + '.nc')
        
    print(fn)
        
    ds = nc.Dataset(fn)
    salt = ds['salt'][0,:,:,:]
    ds.close()
    
    # find the volume and salinity
    
    for seg_name in seg_list:
        net_salt = 0
        full_ji_list = ji_dict[seg_name]
        for ji in full_ji_list:
            net_salt += (salt[:,ji[0],ji[1]] * DV[:,ji[0],ji[1]]).sum()
        volume = v_df.loc[seg_name, 'volume m3']
        mean_salt = net_salt / volume
        # store results
        s_df.loc[dt, seg_name] = mean_salt
        
        if verbose:
            print('%3s: Mean Salinity = %0.4f, Volume = %0.2f km3' %
                (seg_name, mean_salt, volume/1e9))
                
    print('  ** took %0.1f sec' % (time()-tt0))

out_fn = indir + 'daily_segment_salinity.p'
s_df.to_pickle(out_fn)
    
        
    
