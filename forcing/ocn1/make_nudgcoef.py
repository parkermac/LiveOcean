#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:36:47 2017

@author: PM5

Code to make the file of nudging coefficients.
"""

import os
import sys
pth = os.path.abspath('../../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cascadia1', tag='base')
import zrfun

import netCDF4 as nc
import numpy as np

gfn = Ldir['grid'] + 'grid.nc'

G = zrfun.get_basic_info(Ldir['grid'] + 'grid.nc', only_G=True)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

fn = Ldir['grid'] + 'nudgcoef.nc'
# get rid of the old version, if it exists
try:
    os.remove(fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(fn, 'w', format='NETCDF3_CLASSIC')

# create dimensions
foo.createDimension('s_rho', S['N'])
for tag in ['rho', 'u', 'v']:
    foo.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
    foo.createDimension('xi_'+tag, G['lon_'+tag].shape[1])

fld2 = np.zeros_like(G['lon_rho'])
nn = 6
t0 = 1/3
t1 = 1/60

# west 
for i in range(nn):
    for j in range(G['M']):
        tnud = t0 - i*(t0-t1)/nn
        fld2[j, i] = np.max([tnud, fld2[j, i]])
    
# south 
for i in range(G['L']):
    for j in range(nn):
        tnud = t0 - j*(t0-t1)/nn
        fld2[j, i] = np.max([tnud, fld2[j, i]])
    
fld3 = np.zeros((S['N'], G['M'], G['L']))
for i in range(S['N']):
    fld3[i, :, :] = fld2.copy()

# add 2D field data
vv = foo.createVariable('M2_NudgeCoef', float, ('eta_rho', 'xi_rho'))
vv.long_name = '2D momentum inverse nudging coefficients'
vv.units = 'day-1'
vv[:] = fld2

# add 3D field data
vn_dict = {'M3_NudgeCoef': '3D momentum inverse nudging coefficients',
           'tracer_NudgeCoef': 'generic tracer inverse nudging coefficients',
           'temp_NudgeCoef': 'temp inverse nudging coefficients',
           'salt_NudgeCoef': 'salt inverse nudging coefficients'}
for vn in vn_dict.keys():
    vv = foo.createVariable(vn, float, ('s_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = vn_dict[vn]
    vv.units = 'day-1'
    vv[:] = fld3
    
foo.close()
    
"""
From https://www.myroms.org/projects/src/ticket/627

       double M2_NudgeCoef(eta_rho, xi_rho) ;
               M2_NudgeCoef:long_name = "2D momentum inverse nudging coefficients" ;
               M2_NudgeCoef:units = "day-1" ;
               M2_NudgeCoef:coordinates = "xi_rho eta_rho " ;

       double M3_NudgeCoef(s_rho, eta_rho, xi_rho) ;
               M3_NudgeCoef:long_name = "3D momentum inverse nudging coefficients" ;
               M3_NudgeCoef:units = "day-1" ;
               M3_NudgeCoef:coordinates = "xi_rho eta_rho s_rho " ;

       double tracer_NudgeCoef(s_rho, eta_rho, xi_rho) ;
               tracer_NudgeCoef:long_name = "generic tracer inverse nudging coefficients" ;
               tracer_NudgeCoef:units = "day-1" ;
               tracer_NudgeCoef:coordinates = "xi_rho eta_rho s_rho " ;

       double temp_NudgeCoef(s_rho, eta_rho, xi_rho) ;
               temp_NudgeCoef:long_name = "temp inverse nudging coefficients" ;
               temp_NudgeCoef:units = "day-1" ;
               temp_NudgeCoef:coordinates = "xi_rho eta_rho s_rho " ;

       double salt_NudgeCoef(s_rho, eta_rho, xi_rho) ;
               salt_NudgeCoef:long_name = "salt inverse nudging coefficients" ;
               salt_NudgeCoef:units = "day-1" ;
               salt_NudgeCoef:coordinates = "xi_rho eta_rho s_rho " ;
"""

