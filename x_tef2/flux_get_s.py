"""
A tool to extract hourly time series of salinity and volume in the segments.
Now also gets parts of a variance budget.

Performance: takes 1.5 sec per save (3.6 hours per year) on my mac.  This relies
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

import os, sys
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
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2019.07.04')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2019.07.04')
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
outname = Ldir['gtagex'] + '_' + ds0 + '_' + ds1
print(outname +'\n')

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
DX3 = G['DX'].reshape((1,G['M'],G['L']))
DY3 = G['DY'].reshape((1,G['M'],G['L']))

# set input/output location
indir0 = Ldir['LOo'] + 'tef2/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'
#
outdir0 = indir0 + outname + '/'
Lfun.make_dir(outdir0)
outdir = outdir0 + 'flux/'
Lfun.make_dir(outdir)

# load DataFrame of volume and associated dicts
v_df = pd.read_pickle(voldir + 'volumes.p')
bathy_dict = pickle.load(open(voldir + 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(voldir + 'ji_dict.p', 'rb'))

seg_list = list(v_df.index)

testing = False

if testing:
    verbose = True
    seg_list = seg_list[-2:]
else:
    verbose = False

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
s2_df = pd.DataFrame(columns=seg_list)
mix_df = pd.DataFrame(columns=seg_list)
hmix_df = pd.DataFrame(columns=seg_list)
v_df = pd.DataFrame(columns=seg_list)

for fn in fn_list:
    
    tt0 = time()
            
    print(fn)
        
    ds = nc.Dataset(fn)
    salt = ds['salt'][0,:,:,:]
    AKs = ds['AKs'][0,:,:,:]
    KH = float(ds['nl_tnu2'][0].data)
    zeta = ds['zeta'][0,:,:]
    ot = ds['ocean_time'][:]
    ds.close()
    
    # calculate horizontal salinity gradient for hmix
    dsdx = 0*salt
    dsdx[:,:,1:-1] = salt[:,:,2:]-salt[:,:,:-2]
    dsdx = 0.5 * dsdx / DX3[0,:,:]
    
    dsdy = 0*salt
    dsdy[:,1:-1,:] = salt[:,2:,:]-salt[:,:-2,:]
    dsdy = 0.5 * dsdy / DY3[0,:,:]
        
    dt = Lfun.modtime_to_datetime(ot.data[0])
    
    # find the volume and volume-mean salinity
    for seg_name in seg_list:
        
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        z_r, z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S)
        dz = np.diff(z_w, axis=0)
        dzr = np.diff(z_r, axis=0)
        DV = dz * DA3[0,jjj,iii]
        DVR = dzr * DA3[0,jjj,iii]
        volume = DV.sum()
        net_salt = (salt[:,jjj,iii] * DV).sum()
        mean_salt = net_salt/volume
        net_salt2 = (salt[:,jjj,iii] * salt[:,jjj,iii] * DV).sum()
        mean_salt2 = net_salt2/volume
        
        dsdz = (salt[1:,jjj,iii] - salt[:-1,jjj,iii])/dzr
        mix = -2*(AKs[1:-1,jjj,iii] * dsdz * dsdz * DVR).sum()
        
        hmix = -2 * KH * ((dsdx[:,jjj,iii]*dsdx[:,jjj,iii] + dsdy[:,jjj,iii]*dsdy[:,jjj,iii]) * DV).sum()
        
        # store results
        s_df.loc[dt, seg_name] = mean_salt
        s2_df.loc[dt, seg_name] = mean_salt2
        mix_df.loc[dt, seg_name] = mix
        hmix_df.loc[dt, seg_name] = hmix
        v_df.loc[dt, seg_name] = volume
        
        if verbose:
            print('%3s: Mean Salinity = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_salt, volume/1e9))
            print('%3s: Mean Salinity Squared = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_salt2, volume/1e9))
                
    print('  ** took %0.1f sec' % (time()-tt0))
    sys.stdout.flush()

s_out_fn = outdir + 'hourly_segment_salinity.p'
s2_out_fn = outdir + 'hourly_segment_salinity2.p'
mix_out_fn = outdir + 'hourly_segment_mix.p'
hmix_out_fn = outdir + 'hourly_segment_hmix.p'
v_out_fn = outdir + 'hourly_segment_volume.p'

s_df.to_pickle(s_out_fn)
s2_df.to_pickle(s2_out_fn)
mix_df.to_pickle(mix_out_fn)
hmix_df.to_pickle(hmix_out_fn)
v_df.to_pickle(v_out_fn)
    
        
    
