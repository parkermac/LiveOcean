# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

This is designed for NESTING.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing2.sh -g sj0 -t v0 -f ocnN -r backfill -0 20170704 -1 20170704

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

# debugging
testing = True
if testing:
    h_list = h_list[:2]

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
    yes_list = ['zeta', 'ubar', 'vbar'] # only do these 2D fields
    #no_list = ['rho', 'AKv', 'AKs', 'w'] # exclude these 3D fields
    no_list = ['rho', 'AKv', 'AKs', 'w', 'NO3', 'phytoplankton', 'zooplankton',
        'detritus', 'Ldetritus', 'oxygen', 'TIC', 'alkalinity', 'CaCO3']
    for name, v0 in ds0.variables.items():
        
        if (len(v0.dimensions) >= 4 and name not in no_list) or name in yes_list:
        #if name in ['salt']:
            if tt == 0:
                v1 = ds1.createVariable(name, v0.datatype, v0.dimensions)
                v1.time = v0.time
                try:
                    v1.units = v0.units
                except AttributeError:
                    pass # salt has no units
                v1.long_name = v0.long_name
                
            if 'eta_rho' in v0.dimensions:
                tag = 'rho'
            elif 'eta_u' in v0.dimensions:
                tag = 'u'
            elif 'eta_v' in v0.dimensions:
                tag = 'v'
            else:
                print('problem with dimensions')
            X = ds0['lon_' + tag][:]
            Y = ds0['lat_' + tag][:]
            x = G['lon_' + tag]
            y = G['lat_' + tag]
            mask = G['mask_' + tag]
                
            F = ds0[name][:].squeeze()
            if len(F.shape) == 2:
                fx = Ofun.extrap_nearest_to_masked(X, Y, F)
                ff = zfun.interp2(x, y, X, Y, fx)
                fm = np.ma.masked_where(mask==False, ff)
                ds1[name][tt,:,:] = fm
            elif len(F.shape) == 3:
                # For the 3D variables we pull apart the extrapolation and
                # interpolation methods so that the KD Tree and the
                # interpolants are only calculated for the first depth.  This
                # gives a 20x speedup (!) mainly because of the interpolants.
                # This ASSUMES that the mask is the same at all s-levels.
                for iz in range(F.shape[0]):
                    # EXTRAPOLATION
                    fld = F[iz, :, :].squeeze()
                    # do the extrapolation using nearest neighbor
                    fldf = fld.copy() # initialize the "filled" field
                    if iz==0:
                        xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
                        xynew = np.array((X[fld.mask],Y[fld.mask])).T
                        a = cKDTree(xyorig).query(xynew)
                        aa = a[1]
                    fldf[fld.mask] = fld[~fld.mask][aa]
                    fx = fldf.data
                    # INTERPOLATION    
                    if iz==0:
                        # get interpolants
                        xi0, xi1, xf = zfun.get_interpolant(x,X[0,:], extrap_nan=True)
                        yi0, yi1, yf = zfun.get_interpolant(y,Y[:,0], extrap_nan=True)
                    # bi linear interpolation
                    u00 = fx[yi0,xi0]
                    u10 = fx[yi1,xi0]
                    u01 = fx[yi0,xi1]
                    u11 = fx[yi1,xi1]
                    fi = (1-yf)*((1-xf)*u00 + xf*u01) + yf*((1-xf)*u10 + xf*u11)
                    ff = np.reshape(fi, x.shape)
                    fm = np.ma.masked_where(mask==False, ff)
                    ds1[name][tt,iz,:,:] = fm
    print('Hour %d took %0.1f seconds' % (tt, time.time() - tt0))
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


