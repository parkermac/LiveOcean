#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates and a single NetCDF file containing surface height
and other fields, for some time range.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc
import argparse

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

# Command line arguments

# set defaults
gridname = 'cas4'
tag = 'v2'
ex_name = 'lo6biom'
list_type = 'hourly'
# Example of date_string is 2015.09.19
dsf = '%Y.%m.%d'
date_string0 = datetime(2017,1,1).strftime(format=dsf)
date_string1 = datetime(2017,1,3).strftime(format=dsf)

# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str, default=gridname)
parser.add_argument('-t', '--tag', nargs='?', type=str, default=tag)
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default=ex_name)
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default=date_string0)
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default=date_string1)
args = parser.parse_args()
#
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['ex_name'] = args.ex_name
Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']

dt0 = datetime.strptime(args.date_string0, dsf)
dt1 = datetime.strptime(args.date_string1, dsf)
whichyear = dt0.year

# prepare a directory for results
outdir00 = Ldir['parent'] + 'ptools_output/tide/'
Lfun.make_dir(outdir00, clean=False)
outdir0 = outdir00 + 'mod_data/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + Ldir['gtagex'] +'/'
Lfun.make_dir(outdir, clean=False)
# output file
out_name = 'eta_' + str(whichyear) + '.nc'
out_fn = outdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass

if list_type == 'hourly':
    fn_list = []
    dt = dt0
    while dt <= dt1:
        date_string = dt.strftime(format='%Y.%m.%d')
        Ldir['date_string'] = date_string
        f_string = 'f' + Ldir['date_string']
        in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
        if (list_type == 'low_passed') and ('low_passed.nc' in os.listdir(in_dir)):
            fn_list.append(in_dir + 'low_passed.nc') 
        elif list_type == 'hourly':
            for hh in range(2,26):
                hhhh = ('0000' + str(hh))[-4:]
                fn_list.append(in_dir + 'ocean_his_' + hhhh + '.nc')
        dt = dt + timedelta(days=1)
else:
    print('Other list types not implemented yet.')
    
if Ldir['lo_env'] == 'pm_mac':
    pass
    print(fn_list[0])
    print(fn_list[-1])
else:#
    #
    # make some things
    fn = fn_list[0]
    ds = nc.Dataset(fn)
    G = zrfun.get_basic_info(fn, only_G=True)
    h = ds['h'][:]
    ds.close()
    NT = len(fn_list)
    
    ds1 = nc.Dataset(fn_list[0])
    ds2 = nc.Dataset(out_fn, 'w')

    # lists of variables to process
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']

    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        if dname in dlist:
            ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy variables
    for vn in vn_list2:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        vv[:] = ds1[vn][:]
    #
    for vn in ['zeta', 'ocean_time']:
        varin = ds1[vn]
        vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
        vv.long_name = varin.long_name
        vv.units = varin.units
        try:
            vv.time = varin.time
        except AttributeError:
            # ocean_time has no time
            pass

    # copy data
    NT = len(fn_list)
    tt = 0
    for fn in fn_list:
        print(fn)
        sys.stdout.flush()
        ds = nc.Dataset(fn)
        ds2['ocean_time'][tt] = ds['ocean_time'][0].squeeze()
        ds2['zeta'][tt,:,:] = ds['zeta'][0, :, :].squeeze()
        tt += 1
        ds.close()

    ds1.close()
    ds2.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
result_dict['out_fn'] = out_fn
result_dict['date_string0'] = args.date_string0
result_dict['date_string1'] = args.date_string1
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'
print('')
for k in result_dict.keys():
    print('%s: %s' % (k, result_dict[k]))
        





