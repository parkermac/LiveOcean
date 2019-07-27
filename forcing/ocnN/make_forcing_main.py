# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

This is designed for NESTING.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing2.sh -g sj0 -t v0 -f ocnN -r backfill -0 20170704 -1 20170704

and to run from ipython on my mac in LiveOcean/forcing/ocnN/:

run make_forcing_main.py -g sj0 -t v0 -r backfill -d 2017.07.04

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()
import zrfun
import zfun

#%% ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime
import time

#import shutil

import netCDF4 as nc
import numpy as np
from scipy.spatial import cKDTree

opth = os.path.abspath('../ocn4/')
if opth not in sys.path:
    sys.path.append(opth)
import Ofun_nc
import Ofun
from importlib import reload
reload(Ofun_nc)
reload(Ofun)

start_time = datetime.now()
in_dir = Ldir['roms'] + 'output/cas6_v3_lo8b/f' + Ldir['date_string']
nc_dir = Ldir['LOogf_f']

# get grid and S info
G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
# get list of files to work on
h_list_full = os.listdir(in_dir)
h_list = [item for item in h_list_full if 'ocean_his' in item]
h_list.sort()

if Ldir['run_type'] == 'backfill':
    print('Using backfill list')
    h_list = h_list[:25]

# debugging
if Ldir['lo_env'] == 'pm_mac':
    testing = True
else:
    testing = False
    
if testing:
    h_list = h_list[:2]
    
print(h_list[0])
print(h_list[-1])

# get list of times
t_list = []
for h in h_list:
    T = zrfun.get_basic_info(in_dir + '/' + h, only_T=True)
    t_list.append(T['ocean_time'][0])

# name and create output file
clm_fn = nc_dir + 'ocean_clm.nc'
# get rid of the old version, if it exists
try:
    os.remove(clm_fn)
except OSError:
    pass # assume error was because the file did not exist
#ds1 = nc.Dataset(clm_fn, 'w', format='NETCDF3_CLASSIC') # fails
#ds1 = nc.Dataset(clm_fn, 'w', format='NETCDF4_CLASSIC') # works
# the version below works with large files, and it matches the format
# of the history files we are using, so I assume it will work with ROMS.
ds1 = nc.Dataset(clm_fn, 'w', format='NETCDF3_64BIT_OFFSET')

# create dimensions
ds1.createDimension('ocean_time', len(h_list)) # could use None
ds1.createDimension('s_rho', S['N']) # must match ds0.dimensions['N'].size
ds1.createDimension('s_w', S['N'] + 1)
for tag in ['rho', 'u', 'v']:
    ds1.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
    ds1.createDimension('xi_'+tag, G['lon_'+tag].shape[1])    
# add time data    
v1 = ds1.createVariable('ocean_time', float, ('ocean_time',))
v1.units = 'seconds since 1970.01.01 UTC'
v1[:] = np.array(t_list)

# loop over all hours and add processed fields
tt = 0 # hour counter
for h in h_list:
    ds0 = nc.Dataset(in_dir + '/' + h)  
    if S['N'] != ds0.dimensions['N'].size:
        print('Vertical dimensions inconsistent!')
        break
    tt0 = time.time()
    
    for name in ['zeta', 'temp', 'salt']:
        v0 = ds0[name]
        tt00 = time.time()
        if tt == 0:
            v1 = ds1.createVariable(name, v0.datatype, v0.dimensions)
            v1.time = v0.time
            try:
                v1.units = v0.units
            except AttributeError:
                pass # salt has no units
            v1.long_name = v0.long_name
        X = ds0['lon_rho'][:]
        Y = ds0['lat_rho'][:]
        x = G['lon_rho']
        y = G['lat_rho']
        mask = G['mask_rho']
        F = ds0[name][:].squeeze()
        if name == 'zeta':
            fld = F.copy()
            fldf = fld.copy() # initialize the "filled" field
            if tt == 0:
                # fill masked values using nearest neighbor
                xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
                xynew = np.array((X[fld.mask],Y[fld.mask])).T
                a = cKDTree(xyorig).query(xynew)
                aa_rho = a[1]
                # INTERPOLATION    
                # get interpolants
                xi0_rho, xi1_rho, xf_rho = zfun.get_interpolant(x,X[0,:], extrap_nan=True)
                yi0_rho, yi1_rho, yf_rho = zfun.get_interpolant(y,Y[:,0], extrap_nan=True)
            # create the filled field
            fldf[fld.mask] = fld[~fld.mask][aa_rho]
            fx = fldf.data
            # do the bi-linear interpolation
            u00 = fx[yi0_rho,xi0_rho]
            u10 = fx[yi1_rho,xi0_rho]
            u01 = fx[yi0_rho,xi1_rho]
            u11 = fx[yi1_rho,xi1_rho]
            fi = (1-yf_rho)*((1-xf_rho)*u00 + xf_rho*u01) + yf_rho*((1-xf_rho)*u10 + xf_rho*u11)
            ff = np.reshape(fi, x.shape)
            fm = np.ma.masked_where(mask==False, ff)
            ds1[name][tt,:,:] = fm
        elif name in ['temp', 'salt']:
            for iz in range(F.shape[0]):
                # EXTRAPOLATION
                fld = F[iz, :, :].squeeze()
                # do the extrapolation using nearest neighbor
                fldf = fld.copy() # initialize the "filled" field
                fldf[fld.mask] = fld[~fld.mask][aa_rho]
                fx = fldf.data
                # do the bi-linear interpolation
                u00 = fx[yi0_rho,xi0_rho]
                u10 = fx[yi1_rho,xi0_rho]
                u01 = fx[yi0_rho,xi1_rho]
                u11 = fx[yi1_rho,xi1_rho]
                fi = (1-yf_rho)*((1-xf_rho)*u00 + xf_rho*u01) + yf_rho*((1-xf_rho)*u10 + xf_rho*u11)
                ff = np.reshape(fi, x.shape)
                fm = np.ma.masked_where(mask==False, ff)
                ds1[name][tt,iz,:,:] = fm
        print(' -- %s took %0.1f seconds' % (name, time.time() - tt00))
        sys.stdout.flush()

    for name in ['ubar', 'u']:
        v0 = ds0[name]
        tt00 = time.time()
        if tt == 0:
            v1 = ds1.createVariable(name, v0.datatype, v0.dimensions)
            v1.time = v0.time
            v1.units = v0.units
            v1.long_name = v0.long_name
        X = ds0['lon_u'][:]
        Y = ds0['lat_u'][:]
        x = G['lon_u']
        y = G['lat_u']
        mask = G['mask_u']
        F = ds0[name][:].squeeze()
        if name == 'ubar':
            fld = F.copy()
            fldf = fld.copy() # initialize the "filled" field
            if tt == 0:
                # fill masked values using nearest neighbor
                xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
                xynew = np.array((X[fld.mask],Y[fld.mask])).T
                a = cKDTree(xyorig).query(xynew)
                aa_u = a[1]
                # INTERPOLATION    
                # get interpolants
                xi0_u, xi1_u, xf_u = zfun.get_interpolant(x,X[0,:], extrap_nan=True)
                yi0_u, yi1_u, yf_u = zfun.get_interpolant(y,Y[:,0], extrap_nan=True)
            # create the filled field
            fldf[fld.mask] = fld[~fld.mask][aa_u]
            fx = fldf.data
            # do the bi-linear interpolation
            u00 = fx[yi0_u,xi0_u]
            u10 = fx[yi1_u,xi0_u]
            u01 = fx[yi0_u,xi1_u]
            u11 = fx[yi1_u,xi1_u]
            fi = (1-yf_u)*((1-xf_u)*u00 + xf_u*u01) + yf_u*((1-xf_u)*u10 + xf_u*u11)
            ff = np.reshape(fi, x.shape)
            fm = np.ma.masked_where(mask==False, ff)
            ds1[name][tt,:,:] = fm
        elif name == 'u':
            for iz in range(F.shape[0]):
                # EXTRAPOLATION
                fld = F[iz, :, :].squeeze()
                # do the extrapolation using nearest neighbor
                fldf = fld.copy() # initialize the "filled" field
                fldf[fld.mask] = fld[~fld.mask][aa_u]
                fx = fldf.data
                # do the bi-linear interpolation
                u00 = fx[yi0_u,xi0_u]
                u10 = fx[yi1_u,xi0_u]
                u01 = fx[yi0_u,xi1_u]
                u11 = fx[yi1_u,xi1_u]
                fi = (1-yf_u)*((1-xf_u)*u00 + xf_u*u01) + yf_u*((1-xf_u)*u10 + xf_u*u11)
                ff = np.reshape(fi, x.shape)
                fm = np.ma.masked_where(mask==False, ff)
                ds1[name][tt,iz,:,:] = fm
        print(' -- %s took %0.1f seconds' % (name, time.time() - tt00))
        sys.stdout.flush()
        
    for name in ['vbar', 'v']:
        v0 = ds0[name]
        tt00 = time.time()
        if tt == 0:
            v1 = ds1.createVariable(name, v0.datatype, v0.dimensions)
            v1.time = v0.time
            v1.units = v0.units
            v1.long_name = v0.long_name
        X = ds0['lon_v'][:]
        Y = ds0['lat_v'][:]
        x = G['lon_v']
        y = G['lat_v']
        mask = G['mask_v']
        F = ds0[name][:].squeeze()
        if name == 'vbar':
            fld = F.copy()
            fldf = fld.copy() # initialize the "filled" field
            if tt == 0:
                # fill masked values using nearest neighbor
                xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
                xynew = np.array((X[fld.mask],Y[fld.mask])).T
                a = cKDTree(xyorig).query(xynew)
                aa_v = a[1]
                # INTERPOLATION    
                # get interpolants
                xi0_v, xi1_v, xf_v = zfun.get_interpolant(x,X[0,:], extrap_nan=True)
                yi0_v, yi1_v, yf_v = zfun.get_interpolant(y,Y[:,0], extrap_nan=True)
            # create the filled field
            fldf[fld.mask] = fld[~fld.mask][aa_v]
            fx = fldf.data
            # do the bi-linear interpolation
            u00 = fx[yi0_v,xi0_v]
            u10 = fx[yi1_v,xi0_v]
            u01 = fx[yi0_v,xi1_v]
            u11 = fx[yi1_v,xi1_v]
            fi = (1-yf_v)*((1-xf_v)*u00 + xf_v*u01) + yf_v*((1-xf_v)*u10 + xf_v*u11)
            ff = np.reshape(fi, x.shape)
            fm = np.ma.masked_where(mask==False, ff)
            ds1[name][tt,:,:] = fm
        elif name == 'v':
            for iz in range(F.shape[0]):
                # EXTRAPOLATION
                fld = F[iz, :, :].squeeze()
                # do the extrapolation using nearest neighbor
                fldf = fld.copy() # initialize the "filled" field
                fldf[fld.mask] = fld[~fld.mask][aa_v]
                fx = fldf.data
                # do the bi-linear interpolation
                u00 = fx[yi0_v,xi0_v]
                u10 = fx[yi1_v,xi0_v]
                u01 = fx[yi0_v,xi1_v]
                u11 = fx[yi1_v,xi1_v]
                fi = (1-yf_v)*((1-xf_v)*u00 + xf_v*u01) + yf_v*((1-xf_v)*u10 + xf_v*u11)
                ff = np.reshape(fi, x.shape)
                fm = np.ma.masked_where(mask==False, ff)
                ds1[name][tt,iz,:,:] = fm
        print(' -- %s took %0.1f seconds' % (name, time.time() - tt00))
        sys.stdout.flush()
                    
    print('Hour %d took %0.1f seconds' % (tt, time.time() - tt0))
    sys.stdout.flush()
    tt += 1
    ds0.close()
ds1.close()
  
#%% Write other ROMS forcing files
Ofun_nc.make_ini_file(nc_dir)
Ofun_nc.make_bry_file(nc_dir)

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)

ds = nc.Dataset(nc_dir + 'ocean_clm.nc')
ot = ds['ocean_time'][:]
ds.close()
dt0 = Lfun.modtime_to_datetime(ot[0])
dt1 = Lfun.modtime_to_datetime(ot[-1])

result_dict['var_start_time'] = dt0.strftime(time_format)
result_dict['var_end_time'] = dt1.strftime(time_format)

nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
result_dict['result'] = 'success'
for fn in nc_list:
    if os.path.isfile(nc_dir + fn):
        pass
    else:
       result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)


