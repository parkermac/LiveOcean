# -*- coding: utf-8 -*-
"""
Check for min and max of variables in a run.

This is a debugging tool.
"""

#%% setup
import os
import sys
import argparse
from datetime import datetime, timedelta
import netCDF4 as nc4
import numpy as np

sys.path.append(os.path.abspath('../alpha'))
import Lfun

sys.path.append(os.path.abspath('../plotting'))
from importlib import reload
import pfun; reload(pfun)


# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='')
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default='hourly')
args = parser.parse_args()
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# get list of history files to look at

if len(args.date_string0) == 0:
    f_list_raw = os.listdir(Ldir['roms'] + 'output/' + Ldir['gtagex'])
    f_list = [f for f in f_list_raw if f[0]=='f']
    f_list.sort()
    ff = f_list[-1]
    indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + ff + '/'
    fn_list = [indir + 'ocean_his_0001.nc']
else:
    if len(args.date_string1) == 0:
        args.date_string1 = args.date_string0
    fn_list = Lfun.get_fn_list(args.list_type, Ldir, args.date_string0, args.date_string1)
    

print('\n' + 10*'=' + Ldir['gtagex'] + 10*'=')
vn_list = ['u', 'v', 'oxygen', 'TIC', 'alkalinity']
for fn in fn_list:
    print('\n-' + fn.split('/')[-2] + '/' + fn.split('/')[-1])
    ds = nc4.Dataset(fn)
    for vn in vn_list:
        u = ds[vn][0,-1,:,:].squeeze()
        umax, ujmax, uimax, umin, ujmin, uimin = pfun.maxmin(u)
        print('   %10s max=%20.1f (%4d,%4d) min=%20.1f (%4d,%4d)' % (vn, umax, ujmax, uimax, umin, ujmin, uimin))
    ds.close()
    


