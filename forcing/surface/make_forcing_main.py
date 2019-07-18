#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:34:50 2017

@author: PM5

This is the main program for making a SURFACE subset of the daily output.

It creates and a single NetCDF file containing only the surface fields
from all the history files in a given day.

The resulting file for a 3-day cascadia1 forecast with bio and carbon is 327 MB,
vs. 214 MB for EACH of the history files.
This means we are eventually pushing 47x less to azure.

For testing on my mac run in ipython as
run make_forcing_main.py -d 2017.05.18

2019.07.18 Shortened the list of variables to only those served by EDS Viewer.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime
start_time = datetime.now()
import netCDF4 as nc
import zrfun
import numpy as np

print(' - Creating surface  file for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']

#%% create the surface NetCDF file

# input files
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = []
for item in fn_list_raw:
    if 'ocean_his' in item and '.nc' in item:
        fn_list.append(in_dir + item)
fn_list.sort()

# shorten the list to only be every 2 hours
fn_list = fn_list[::2]

#%% make z
# fn = fn_list[0]
# ds = nc.Dataset(fn)
# S = zrfun.get_basic_info(fn, only_S=True)
# h = ds['h'][:]
# z = zrfun.get_z(h, 0*h, S, only_rho=True)
# z0 = z[-1,:,:].squeeze()
# ds.close()

#%% Initialize the multi-file input dataset
ds1 = nc.MFDataset(fn_list)

#%% make surface velocity
u0 = ds1['u'][:, -1, :, :].squeeze()
v0 = ds1['v'][:, -1, :, :].squeeze()
u = np.nan * ds1['salt'][:, -1, :, :].squeeze()
v = u.copy()
u[:, :, 1:-1] = (u0[:, :, 1:] + u0[:, :, :-1])/2
v[:, 1:-1, :] = (v0[:, 1:, :] + v0[:, :-1, :])/2

# output files
out_name = 'ocean_surface.nc'
out_fn = in_dir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list2t = ['Uwind', 'Vwind', 'ocean_time']
vn_list3t = ['salt', 'temp']

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
for vn in vn_list2t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.long_name = varin.long_name
    vv.units = varin.units
    try:
        vv.time = varin.time
    except AttributeError:
        # ocean_time has no time
        pass
    vv[:] = ds1[vn][:]
#
for vn in vn_list3t:
    do_var = True
    try:
        varin = ds1[vn]
    except IndexError:
        # designed so that it still works when we don't have a variable from this list
        # e.g. when there is no bio or carbon
        do_var = False
        print(' - Variable not found: ' + vn)
    if do_var==True:
        dd = tuple([d for d in varin.dimensions if d != 's_rho'])
        vv = ds2.createVariable(vn, varin.dtype, dd)
        if vn == 'PH':
            vv.long_name = 'pH'
        elif vn == 'ARAG':
            vv.long_name = 'Aragonite Saturation State'
        else:
            vv.long_name = varin.long_name
        try:
            vv.units = varin.units
        except AttributeError:
            # salt has no units
            pass
        vv.time = varin.time
        vv[:] = ds1[vn][:, -1, :, :].squeeze()

# Add derived variables
# vv = ds2.createVariable('z', float, ('eta_rho', 'xi_rho'))
# vv.long_name = 'z position closest to free surface for 3D variables '
# vv.units = 'meter'
# vv[:] = z0
#
vv = ds2.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'eastward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = u
#
vv = ds2.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'northward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = v

ds1.close()
ds2.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
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

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)

