"""
Code to test getting a single time slice of HYCOM FMRC fields.
"""

import os, sys
sys.path.append(os.path.abspath('../../alpha/'))
import Lfun
Ldir = Lfun.Lstart()

import netCDF4 as nc
from urllib.request import urlretrieve
from datetime import datetime, timedelta
import numpy as np

dt = datetime(2020,3,14)
dstr = dt.strftime('%Y-%m-%d-T00:00:00Z')

outdir = Ldir['LOo'] + 'test/'
Lfun.make_dir(outdir)
fn_out = outdir + 'fmrc_oneday_test.nc'

try:
    os.remove(fn_out)
except OSError:
    pass # assume error was because the file did not exist

url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'
    + '?var=surf_el&north=53&south=39&west=229&east=239'
    + '&time_start='+dstr+'&time_end='+dstr
    + '&addLatLon=true&accept=netcdf4')

(a,b) = urlretrieve(url,fn_out)

ds = nc.Dataset(fn_out)

# get time info for the forecast
t = ds['time'][0]
if isinstance(t, np.ma.MaskedArray):
    th = t.data
else:
    th = t
tu = ds['time'].units
# e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
# Warning: Brittle code below!
ymd = tu.split()[2]
hmss = tu.split()[3]
hms = hmss.split('.')[0]
hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
this_dt = hycom_dt0 + timedelta(days=(th/24))
print('Target time = ' + dstr)
print('Actual time = ' + this_dt.strftime('%Y-%m-%d-T00:00:00Z'))
