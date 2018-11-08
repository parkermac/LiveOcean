"""
This is the main program for making the UBC output file.

For testing on my mac run in ipython as
run make_forcing_main.py -d 2017.05.18

Runs in under a second!
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

# modified 2018.11.08 by PM to use a lat-lon range instead of indices, prompted
# by the change from cascadia1_base_lobio5 to cas4_v2_lo6biom.

# Also makes a low-passed version of the files, instead of relying on the existence
# of the low-passed file.

testing=False

pth = os.path.abspath('../../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import Lfun
import zrfun

from importlib import reload
import UBC_subdomain
reload(UBC_subdomain)

import numpy as np
import shutil

# generate a list of all the files
in_dir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
    '/f' + Ldir['date_string'] + '/')
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir + ff) for ff in fn_list_raw if 'ocean_his' in ff]
fn_list.sort()
if len(fn_list) != 73:
    print('Warning: We have %d history files' % (len(fn_list)))
    
if testing == True:
    fn_list = fn_list[:7]

G = zrfun.get_basic_info(fn_list[0], only_G=True)

lon_vec = G['lon_rho'][0,:]
lat_vec = G['lat_rho'][:,0]

# test values for cascadia1
# XBS = [55, 80]  # x-limits
# YBS = [295, 325]  # y-limits

# Bounds of subdomain (rho grid) from Susan Allen 2018_11
lon0 = -125.016452048434
lon1 = -124.494612925929
lat0 = 48.3169463809796
lat1 = 48.7515055163539

# find indices containing these bonds
i0, i1, ifr = zfun.get_interpolant(np.array([lon0,lon1]), lon_vec, extrap_nan=True)
j0, j1, jfr = zfun.get_interpolant(np.array([lat0,lat1]), lat_vec, extrap_nan=True)
if np.nan in ifr:
    print('Warning: nan in ifr')
if np.nan in jfr:
    print('Warning: nan in jfr')
XBS = [i0[0], i1[1]]
YBS = [j0[0], j1[1]]
print(XBS)
print(YBS)

from datetime import datetime, timedelta
start_time = datetime.now()

# the output file name
out_fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
    '/f' + Ldir['date_string'] + '/low_passed_UBC.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
    
# make a temporary place for the output
temp_dir = in_dir + 'ubc_temp/'
Lfun.make_dir(temp_dir, clean=True)

# do the extraction for the whole list
UBC_subdomain.get_UBC_subdomain(fn_list, temp_dir, XBS, YBS)

# prepare for the low_pass
out_fn_list_raw = os.listdir(temp_dir)
out_fn_list = [(temp_dir + ff) for ff in out_fn_list_raw if 'ocean_his' in ff]
out_fn_list.sort()
# create the filter
out_fn_list_short = out_fn_list[1:-1] # trim to get to 71 for a forecast
nf = len(out_fn_list_short)
if nf == 71:
    print(' - Using Godin filter')
    filt0 = zfun.godin_shape()
else:
    print(' - Using Hanning filter for list length = ' + str(nf))
    filt0 = zfun.hanning_shape(nf)
# and do the low pass
zrfun.roms_low_pass(out_fn_list_short, out_fn, filt0, exclude=[])

# get rid of the temp directory
shutil.rmtree(temp_dir, ignore_errors=True)

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

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)