#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Extract as-run river time series.

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


# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.01.01')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2018.12.31')
args = parser.parse_args()

# save some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)

dt0 = datetime.strptime(args.date_string0, '%Y.%m.%d')
dt1 = datetime.strptime(args.date_string1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1

mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, '%Y.%m.%d'))
    mdt = mdt + timedelta(days=1)


# get some initial info
mds = mds_list[0]
fn = Ldir['LOo'] + Ldir['gtag'] + '/f' + mds + '/riv2/rivers.nc'
ds = nc.Dataset(fn)
rn = ds['river_name'][:]
NR = rn.shape[1]
riv_name_list = []
for ii in range(NR):
    a = rn[:,ii]
    r = []
    for l in a:
        r.append(l.decode())
    rr = ''.join(r)
    #print('%d %s' % (ii, rr))
    riv_name_list.append(rr)
ds.close()

NT = len(mds_list)

Qr = np.nan * np.ones((NT, NR))
tt = 0
for mds in mds_list:
    fn = Ldir['LOo'] + Ldir['gtag'] + '/f' + mds + '/riv2/rivers.nc'
    ds = nc.Dataset(fn)
    # the river transport is given at noon of a number of days surrounding the forcing date
    # here we find the index of the time for today
    RT = ds['river_time'][:]
    ii = 0
    for rt in RT:
        rdt = Lfun.modtime_to_datetime(rt)
        rds = datetime.strftime(rdt, '%Y.%m.%d')
        if rds == mds:
            #print('Match at ii=%d, %s' % (ii, rds))
            Qr[tt,:] = ds['river_transport'][ii,:]
        ii += 1
    ds.close()
    tt += 1
    
Qr = np.abs(Qr)

# now interpolate to make the time and Qr at midnight of each day
Qind = pd.date_range(dt0, dt1+timedelta(days=1),freq='D')
Qr = np.concatenate((np.reshape(Qr[0,:], (1,NR)),Qr,np.reshape(Qr[-1,:], (1,NR))), axis=0)
Qr = Qr[:-1,:] + np.diff(Qr, axis=0)/2

df = pd.DataFrame(index=Qind, columns=riv_name_list, data=Qr)

outdir = Ldir['LOo'] + 'river/'
Lfun.make_dir(outdir)
fn = outdir + Ldir['gtag'] + '_' + args.date_string0 + '_' + args.date_string1 + '.p'

df.to_pickle(fn)
    
