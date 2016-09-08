# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 15:44:54 2016

@author: PM5

This is the main program for making the TIDE forcing file.

It is designed for Analytical runs.

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()


#%% ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime
import netCDF4 as nc
import zrfun
import numpy as np

start_time = datetime.now()

out_fn = Ldir['LOogf_f'] + 'tides.nc'
grid_fn = Ldir['grids'] + Ldir['gridname'] + '/grid.nc'

[G] = zrfun.get_basic_info(grid_fn, getS=False, getT=False)

dst = nc.Dataset(out_fn, 'w', format='NETCDF3_CLASSIC')

# get sizes and initialize the NetCDF file
# specify the constituents to use
cons_list =  ['m2', 's2']
ncons = len(cons_list)

lon_rho = G['lon_rho']
lat_rho = G['lat_rho']
mask_rho = G['mask_rho'];
ny, nx = lon_rho.shape

dst.createDimension('tide_period', ncons)
dst.createDimension('slen', 10)
dst.createDimension('xi_rho', nx)
dst.createDimension('eta_rho', ny)

v_var = dst.createVariable('constituent_name', 'c', ('slen', 'tide_period'))
rr = 0
for cons in cons_list:
    print('- ' + cons)
    cc = 0
    for ch in cons:
        v_var[cc, rr] = ch
        cc += 1
    rr += 1

v_var = dst.createVariable('tide_period', float, ('tide_period'))
v_var[:] = [12.42, 12.0]
v_var.units = 'hours'

v_var = dst.createVariable('tide_Eamp', float, ('tide_period', 'eta_rho', 'xi_rho'))
v_var[0, :, :] = np.ones_like(lon_rho).astype(float)
v_var[1, :, :] = np.zeros_like(lon_rho).astype(float)

v_var = dst.createVariable('tide_Ephase', float, ('tide_period', 'eta_rho', 'xi_rho'))
v_var[:] = 0.

dst.close()

"""
    % and make final adjustments before writing to arrays
    tide_period(ii) = 2*pi/(3600*om); % hours
    tide_Eamp(ii,:,:) = pf*Eamp; % m
    tide_Ephase(ii,:,:) = Ephase - 180*ph/pi - 180*pu/pi; % deg
    tide_Cangle(ii,:,:) = Cangle; % deg
    tide_Cphase(ii,:,:) = Cphase - 180*ph/pi - 180*pu/pi; % deg
    tide_Cmax(ii,:,:) = pf*Cmax/100; % m s-1
    tide_Cmin(ii,:,:) = pf*Cmin/100; % m s-1
"""

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

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
