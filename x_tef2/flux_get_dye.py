"""
A tool to extract hourly time series of dye and volume in the segments.

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
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8dye')
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
DXu = (G['DX'][:,1:]+G['DX'][:,:-1])/2
DX3u = DXu.reshape((1,G['M'],G['L']-1))
DYv = (G['DY'][1:,:]+G['DY'][:-1,:])/2
DY3v = DYv.reshape((1,G['M']-1,G['L']))

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

dye_df = pd.DataFrame(columns=seg_list)
v_df = pd.DataFrame(columns=seg_list)

for fn in fn_list:
    
    tt0 = time()
            
    print(fn)
        
    ds = nc.Dataset(fn)
    dye = ds['dye_01'][0,:,:,:]
    zeta = ds['zeta'][0,:,:]
    ot = ds['ocean_time'][:]
    ds.close()
            
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
        net_dye = (dye[:,jjj,iii] * DV).sum()
        mean_dye = net_dye/volume
        
        
        # store results
        dye_df.loc[dt, seg_name] = mean_dye
        v_df.loc[dt, seg_name] = volume
        
        if verbose:
            print('%3s: Mean Dye = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_dye, volume/1e9))
                
    print('  ** took %0.1f sec' % (time()-tt0))
    sys.stdout.flush()

dye_out_fn = outdir + 'hourly_segment_dye.p'
dye_v_out_fn = outdir + 'hourly_segment_dye_volume.p'

dye_df.to_pickle(dye_out_fn)
v_df.to_pickle(dye_v_out_fn)
    
        
    
