#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This extracts selected fields from the model forcing.  Initially coded
to look at the SSH in the bry file, as part of the NSoG problem.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc
import argparse
import pandas as pd

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

# Command line arguments

# defaults
list_type = 'hourly' # hourly, daily, low_passed
layer_name = 'zeta' # name to control choices in code

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.01.01')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.01.02')
args = parser.parse_args()

# get time limits
ds0 = args.date_string0; ds1 = args.date_string1
dt0 = datetime.strptime(ds0, '%Y.%m.%d'); dt1 = datetime.strptime(ds1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1

Ldir = Lfun.Lstart(args.gridname, args.tag)

# get some grid info
grid_fn = Ldir['grid'] + 'grid.nc'
dsg = nc.Dataset(grid_fn)
mask_north = dsg['mask_rho'][-1,:].squeeze()
dsg.close()

# make sure the output directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = outdir00 + 'misc/'
Lfun.make_dir(outdir0)
outdir = (outdir0 + 'bry_' + Ldir['gtag'] + '_' + ds0 + '_' + ds1 + '/')
Lfun.make_dir(outdir)

date_list = Lfun.date_list_utility(dt0, dt1)

# name output file
out_fn = (outdir + 'zeta_north.p')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass

# decide which ocn to use
if Ldir['gtag'] == 'cas4_v2':
    ocn = 'ocn2'
elif Ldir['gtag'] == 'cas5_v3':
    ocn = 'ocn4'
    
# initialize output vectors
NT = len(date_list)
zeta_north_ser = pd.Series()
for date_str in date_list:
    in_fn = Ldir['LOo'] + Ldir['gtag'] + '/f' + date_str + '/' + ocn + '/ocean_bry.nc'
    try:
        ds = nc.Dataset(in_fn)
        zeta_north = ds['zeta_north'][0,:].squeeze()
        zeta_north_mean = zeta_north[mask_north == 1].mean()
        dt = datetime.strptime(date_str, '%Y.%m.%d')
        zeta_north_ser[dt] = zeta_north_mean
        ds.close()
    except:
        pass
        
# save output
zeta_north_ser.to_pickle(out_fn)
