""""
Code to test getting hycom files using the new FMRC_best file.
"""

import os
import netCDF4 as nc

from urllib.request import urlretrieve
from urllib.error import URLError
from socket import timeout
import time
from datetime import datetime, timedelta

testing = True

# specify output file
fn_out = 'test.nc'
# get rid of the old version, if it exists
try:
    os.remove(fn_out)
except OSError:
    pass # assume error was because the file did not exist

# specify time limits
# get today's  date string in LO format
dstr00 = datetime.now().strftime('%Y.%m.%d')
print('Working on ' + dstr00)
# get time limits for forecast
dt00 = datetime.strptime(dstr00, '%Y.%m.%d')
dt0 = dt00 - timedelta(days=2)
if testing == True:
    dt1 = dt00 - timedelta(days=1)
else:
    dt1 = dt00 + timedelta(days=5)
# put them in ncss format
dstr0 = dt0.strftime('%Y-%m-%d-T00:00:00Z')
dstr1 = dt1.strftime('%Y-%m-%d-T00:00:00Z')
print('- dt0 = ' + dstr0)
print('- dt1 = ' + dstr1)

# specify spatial limits
#aa = [-129, -121, 39, 51]
# # find indices of a sub region
# aa = hfun.get_extraction_limits()
north = 51
south = 39
west = 231 #-129 + 360
east = 238 #-122 + 360

if testing == True:
    var_list = 'surf_el'
    #var_list = 'surf_el,salinity'
else:
    var_list = 'surf_el,water_temp,salinity,water_u,water_v'

# create the request url
url = ('http://ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0/data/forecasts/FMRC_best.ncd'+
    '?var='+var_list +
    '&north='+str(north)+'&south='+str(south)+'&west='+str(west)+'&east='+str(east) +
    '&time_start='+dstr0+'&time_end='+dstr1+'&timeStride=8' +
    '&addLatLon=true&accept=netcdf4')

# get the data and save as a netcdf file
counter = 1
got_file = False
while (counter <= 3) and (got_file == False):
    print('Attempting to get data, counter = ' + str(counter))
    tt0 = time.time()
    try:
        (a,b) = urlretrieve(url,fn_out)
        # a is the output file name
        # b is a message you can see with b.as_string()
    except URLError as e:
        if hasattr(e, 'reason'):
            print(' *We failed to reach a server.')
            print(' -Reason: ', e.reason)
        elif hasattr(e, 'code'):
            print(' *The server couldn\'t fulfill the request.')
            print(' -Error code: ', e.code)
    except timeout:
        print(' *Socket timed out')
    else:
        got_file = True
        print(' Worked fine')
    print(' -took %0.1f seconds' % (time.time() - tt0))
    counter += 1

# check results
ds = nc.Dataset('test.nc')
print('\nVariables:')
for vn in ds.variables:
    print('- '+vn)
print('Times:')
htime = ds['time'][:]
# get time info for the forecast
t = ds['time'][:]
tu = ds['time'].units
# e.g. 'hours since 2018-11-20 12:00:00.000 UTC'
# Warning: Brittle code below!
ymd = tu.split()[2]
hmss = tu.split()[3]
hms = hmss.split('.')[0]
hycom_dt0 = datetime.strptime(ymd + ' ' + hms, '%Y-%m-%d %H:%M:%S')
# make a list of datetimes that are in the forecast
hycom_dt_list = [] # datetimes
hycom_iit_list = [] # indices
iit = 0
for th in t: # assume t is sorted
    this_dt = hycom_dt0 + timedelta(th/24)
    hycom_dt_list.append(this_dt)
    print('- '+str(this_dt))
    hycom_iit_list.append(iit)
    iit += 1
ds.close()
