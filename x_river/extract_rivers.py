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
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2018.09.01')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2018.09.03')
args = parser.parse_args()

# save some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)

fn = Ldir['LOo'] + Ldir['gtag'] + '/f' + args.date_string0 + '/riv2/rivers.nc'

ds = nc.Dataset(fn)

# for vn in ds.variables:
#     print(vn)

    
rn = ds['river_name'][:]
NR = rn.shape[1]
for ii in range(NR):
    a = rn[:,ii]
    r = []
    for l in a:
        r.append(l.decode())
    rr = ''.join(r)
    print('%d %s' % (ii, rr))
