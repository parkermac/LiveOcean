# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 08:24:42 2016

@author: PM5

This is the main program for making the OCN forcing file.

It is for Analytical runs.

*******************************

To test in python on mac:

cd ~/Documents/LiveOcean/forcing/ocnA

run make_forcing_main.py -g aestus1 -t A1 -r backfill -d 2013.01.01

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

#%% ****************** CASE-SPECIFIC CODE *****************

import zfun
import zrfun
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
out_dir = Ldir['LOogf_f']

#%% make climatology

# get grid and S info
G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

# name output file
clm_fn = out_dir + 'ocean_clm.nc'
# get rid of the old version, if it exists
try:
    os.remove(clm_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(clm_fn, 'w', format='NETCDF3_CLASSIC')

# create dimensions
for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
    foo.createDimension(vn+'_time', 2)
foo.createDimension('s_rho', S['N'])
for tag in ['rho', 'u', 'v']:
    foo.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
    foo.createDimension('xi_'+tag, G['lon_'+tag].shape[1])

# add time data
dtf = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
dt0 = dtf
dt1 = dtf + timedelta(days=5)
dt0m = Lfun.datetime_to_modtime(dt0)
dt1m = Lfun.datetime_to_modtime(dt1)
for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
    vv = foo.createVariable(vn+'_time', float, (vn+'_time',))
    vv.units = 'seconds since 1970.01.01 UTC'
    vv[:] = np.array([dt0m, dt1m])

# add 2d field data
vv = foo.createVariable('zeta', float, ('zeta_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'sea surface height climatology'
vv.units = 'meter'
vv[:] = 0
vv = foo.createVariable('ubar', float, ('v2d_time', 'eta_u', 'xi_u'))
vv.long_name = 'vertically averaged u-momentum climatology'
vv.units = 'meter second-1'
vv[:] = 0
vv = foo.createVariable('vbar', float, ('v2d_time', 'eta_v', 'xi_v'))
vv.long_name = 'vertically averaged v-momentum climatology'
vv.units = 'meter second-1'
vv[:] = 0

# add 3d field data
vv = foo.createVariable('u', float, ('v3d_time', 's_rho', 'eta_u', 'xi_u'))
vv.long_name = 'u-momentum component climatology'
vv.units = 'meter second-1'
vv[:] = 0
vv = foo.createVariable('v', float, ('v3d_time', 's_rho', 'eta_v', 'xi_v'))
vv.long_name = 'v-momentum component climatology'
vv.units = 'meter second-1'
vv[:] = 0
vv = foo.createVariable('salt', float, ('salt_time', 's_rho', 'eta_rho', 'xi_rho'))
vv.long_name = 'salinity climatology'
vv.units = 'PSU'
vv[:] = 35
vv[:, :, :, 75:] = 0
vv = foo.createVariable('temp', float, ('temp_time', 's_rho', 'eta_rho', 'xi_rho'))
vv.long_name = 'potential temperature climatology'
vv.units = 'Celsius'
vv[:] = 10

foo.close()
#print('============= clim ===========================================')
#zfun.ncd(clm_fn) # testing

#%% Initial condition, copied from first time of clm_fn

ds1 = nc.Dataset(clm_fn, mode='r')
# name output file
ini_fn = out_dir + 'ocean_ini.nc'
# get rid of the old version, if it exists
try:
    os.remove(ini_fn)
except OSError:
    pass # assume error was because the file did not exist
ds2 = nc.Dataset(ini_fn, 'w', format='NETCDF3_CLASSIC')

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    if 'time' in dname:
        ds2.createDimension(dname, 1)
    else:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy variables
for v_name, varin in ds1.variables.items():
    outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
    # Copy variable attributes, {} is a dict comprehension, cool!
    outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
    if varin.ndim > 1:
        outVar[:] = varin[0,:]
    else:
        outVar[:] = varin[0]

ds1.close()
ds2.close()
#print('============= ini ===========================================')
#zfun.ncd(ini_fn) # testing

#%% Boundary conditions, copied from edges of clm_fn

ds1 = nc.Dataset(clm_fn, mode='r')
# name output file
bry_fn = out_dir + 'ocean_bry.nc'
# get rid of the old version, if it exists
try:
    os.remove(bry_fn)
except OSError:
    pass # assume error was because the file did not exist
ds2 = nc.Dataset(bry_fn, 'w', format='NETCDF3_CLASSIC')

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy parts of variables
for v_name, varin in ds1.variables.items():

    if varin.ndim in [3,4]: 
        for bname in ['north', 'south', 'east', 'west']:
            outname = v_name + '_' + bname
            if bname in ['north', 'south']:
                outdims = tuple([item for item in varin.dimensions if item[:3] != 'eta'])
            elif bname in ['west', 'east']:
                outdims = tuple([item for item in varin.dimensions if item[:2] != 'xi'])
            outVar = ds2.createVariable(outname, varin.datatype, outdims)    
            outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})            
            print(outname)
            print(outdims)
            print(varin.ndim)
            print(varin.shape)
            if varin.ndim == 4:
                if bname == 'north':
                    outVar[:] = varin[:,:,-1,:]
                    print('hi4')
                elif bname == 'south':
                    outVar[:] = varin[:,:,0,:]
                elif bname == 'east':
                    outVar[:] = varin[:,:,:,-1]
                elif bname == 'west':
                    outVar[:] = varin[:,:,:,0]
            elif varin.ndim == 3:
                if bname == 'north':
                    print('hi3')
                    outVar[:] = varin[:,-1,:]
                elif bname == 'south':
                    outVar[:] = varin[:,0,:]
                elif bname == 'east':
                    outVar[:] = varin[:,:,-1]
                elif bname == 'west':
                    outVar[:] = varin[:,:,0]
    else:
        outname = v_name
        outdims = tuple([item for item in varin.dimensions])
        outVar = ds2.createVariable(outname, varin.datatype, outdims)    
        outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
        outVar[:] = varin[:]

    
    
ds1.close()
ds2.close()
print('============= bry ===========================================')
zfun.ncd(bry_fn) # testing

# ************** END CASE-SPECIFIC CODE *****************

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
result_dict['var_start_time'] = dt0.strftime(time_format)
result_dict['var_end_time'] = dt1.strftime(time_format)
if os.path.isfile(bry_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)


