"""
Code to test getting and processing the HYOCM fields using ncks.

See email thread with Michael McDonald around 2020.06.23 in LiveOcean/HYCOM.

At the linux command line use "ncks -help" to see usage notes.

RESULT: appears to work reliably, and takes about 4.5 mintues on my mac.
"""

import os, sys
sys.path.append(os.path.abspath('../../alpha/'))
import Lfun
Ldir = Lfun.Lstart()

import netCDF4 as nc
from datetime import datetime, timedelta
import numpy as np
from time import time
import subprocess

dt_now = datetime.now()
dt00 = datetime(dt_now.year,dt_now.month,dt_now.day) # strip HMS

testing = False

# set time range
if testing:
    # shorter
    dt0 = dt00
    dt1 = dt0 + timedelta(days=1)
else:
    # full forecast
    dt0 = dt00 - timedelta(days=2)
    dt1 = dt0 + timedelta(days=7)
dstr0 = dt0.strftime('%Y-%m-%dT00:00') 
dstr1 = dt1.strftime('%Y-%m-%dT00:00')

# output info
outdir = Ldir['LOo'] + 'test/'
Lfun.make_dir(outdir)
fn = outdir + 'ncks_test.nc'
# try:
#     os.remove(fn)
# except OSError:
#     pass
    
# use subprocess.call() to execute the ncks command
if testing:
    vstr = 'surf_el,depth'
else:
    vstr = 'surf_el,water_temp,salinity,water_u,water_v,depth'
cmd_list = ['ncks',
    #'-D8', # debug level
    '-d', 'time,'+dstr0+','+dstr1+',8',
    '-d', 'lon,229.,239.,1',
    '-d','lat,39.,53.,1',
    '-v',vstr,
    'https://tds.hycom.org/thredds/dodsC/GLBy0.08/latest',
    '-4', # NetCDF4 output
    '-O', fn]  # -O means overwrite output

# run ncks
tt0 = time()
ret = subprocess.call(cmd_list)
print('Time to get file = %0.2f sec' % (time()-tt0))
print('Return code = ' + str(ret) + ' (0=success)')

def get_hdt_list(fn):
    # function to get the times of a HYCOM ncks extraction as a list of datetimes
    ds = nc.Dataset(fn)
    # get info for the forecast
    t = ds['time'][:]
    if isinstance(t, np.ma.MaskedArray):
        th_vec = t.data
    else:
        th_vec = t
    tu = ds['time'].units
    # e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
    ymd = tu.split()[2]
    hmss = tu.split()[3]
    hms = hmss.split('.')[0]
    hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
    hdt_list = []
    for th in th_vec:
        this_dt = hycom_dt0 + timedelta(days=(th/24))
        hdt_list.append(this_dt)
    ds.close()
    return hdt_list
    
# check output
# - times
hdt_list = get_hdt_list(fn)
print('\nTarget time range = ' + dstr0 + ' to ' + dstr1)
for hdt in hdt_list:
    print(' - Actual time = ' + hdt.strftime('%Y-%m-%d-T00:00:00Z'))
print('')
# - variables
ds = nc.Dataset(fn)
for vn in ds.variables:
    print(vn + ' [' + ds[vn].units+ '] ' + str(ds[vn].shape))
ds.close()

# next split it into individual files as expected by the later processing
print('')
print('Split up into separate output files:')
NT = len(hdt_list)
for ii in range(NT):
    hdt = hdt_list[ii]
    iis = str(ii)
    fn_out = outdir + 'h'+ hdt.strftime('%Y.%m.%d') + '.nc'
    cmd_list = ['ncks',
        '-d', 'time,'+iis+','+iis,
        '-O', fn, fn_out]
    ret = subprocess.call(cmd_list)
    this_hdt = get_hdt_list(fn_out)[0]
    print(fn_out + ': actual time = ' + str(this_hdt))
    
print('\nVaribles of last split file:')
# check variables of the last split file
ds = nc.Dataset(fn_out)
for vn in ds.variables:
    print(vn + ' [' + ds[vn].units+ '] ' + str(ds[vn].shape))
ds.close()

