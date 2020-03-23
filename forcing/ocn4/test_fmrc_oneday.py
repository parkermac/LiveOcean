"""
Code to test getting a single time slice of HYCOM FMRC fields.
"""

import os, sys
sys.path.append(os.path.abspath('../../alpha/'))
import Lfun
Ldir = Lfun.Lstart()

import netCDF4 as nc
from urllib.request import urlretrieve
import requests
from datetime import datetime, timedelta
import numpy as np
from time import time

dt = datetime.now() + timedelta(days=1)
dstr = dt.strftime('%Y-%m-%d-T00:00:00Z')

outdir = Ldir['LOo'] + 'test/'
Lfun.make_dir(outdir)
fn_out = outdir + 'fmrc_oneday_test.nc'
fn_out2 = outdir + 'fmrc_oneday_test2.nc'

try:
    os.remove(fn_out)
except OSError:
    pass # assume error was because the file did not exist

try:
    os.remove(fn_out2)
except OSError:
    pass # assume error was because the file did not exist

if False:
    # just get SSH
    url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'
        + '?var=surf_el'
        + '&north=53&south=39&west=229&east=239'
        + '&time'+dstr
        + '&addLatLon=true&accept=netcdf4')
else:
    # get all variables
    url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'
        + '?var=surf_el,water_temp,salinity,water_u,water_v'
        + '&north=53&south=39&west=229&east=239'
        + '&time=' + dstr
        + '&addLatLon=true&accept=netcdf4')

def get_time(fn):
    ds = nc.Dataset(fn)
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
    ds.close()

tt0 = time()
(a,b) = urlretrieve(url,fn_out)
print('\nurlretrieve took %0.1f seconds' % (time()-tt0))
get_time(fn_out)

tt0 = time()
if False:
    # streaming, chunked version
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(fn_out2, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
else:
    # regular version
    r = requests.get(url)
    with open(fn_out2,'wb') as f:
        f.write(r.content)
print('\nrequests.get took %0.1f seconds' % (time()-tt0))
get_time(fn_out2)
