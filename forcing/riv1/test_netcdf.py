# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 13:47:36 2016

@author: PM5

Code to experiment with writing NetCDF files.

"""

import netCDF4 as nc

fn = '/Users/PM5/Desktop/test.nc'

foo = nc.Dataset(fn, 'w', format='NETCDF3_CLASSIC')

rn_list = ['skagit', 'columbia']

N = len(rn_list)
SL = 50

foo.createDimension('river', N)
foo.createDimension('slen', SL)

v_var = foo.createVariable('river_name', 'c', ('river', 'slen'))

print('\nInput:')
rr = 0
for rn in rn_list:
    print('- ' + rn)
    cc = 0
    for ch in rn:
        v_var[rr,cc] = ch
        cc += 1
    rr += 1

foo.close()

#%% check results

ds = nc.Dataset(fn)

river_names = ds['river_name'][:]

print('\nOutput:')
for rn_arr in river_names:
    RN = ''
    for ch in rn_arr:
        # we use .decode() to get rid of the b prefix on A
        RN = RN + ch.decode()
    print('- ' + RN)

ds.close()
