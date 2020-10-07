"""
This is the main program for making the RIV forcing file.

It it just for adding dye values, and copies form existing files.

to test from mac in ipython

run make_forcing_main.py -g cas6 -t v3 -r backfill -d 2019.07.04

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)

from importlib import reload
import forcing_functions as ffun
reload(ffun)
Ldir, Lfun = ffun.intro()


# ****************** CASE-SPECIFIC CODE *****************
import netCDF4 as nc
import shutil
from datetime import datetime
start_time = datetime.now()

# copy rivers.nc file from riv2 to riv2dye
out_frc = Ldir['frc']
in_frc = out_frc.replace('dye','')
in_fn = Ldir['LOogf'] + in_frc + '/rivers.nc'
out_fn = Ldir['LOogf_f'] + 'rivers.nc'
print('in_fn = ' + in_fn)
print('out_fn = ' + out_fn)
if os.path.exists(out_fn):
    os.remove(out_fn)
shutil.copyfile(in_fn, out_fn)

# edit that file
ds = nc.Dataset(out_fn, 'a')
vn = 'river_dye_01'
if vn not in ds.variables:
    vv = ds.createVariable(vn, float, ('river_time', 's_rho', 'river'))
    vv.long_name = 'river dye'
    vv.time = "river_time"
    vv.units = "kg m-3"
    vv[:] = 0 * ds['river_salt'][:]
else:
    print('** WARNING ' + vn + ' is already in the output file! **')
ds.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
# result_dict['var_start_time'] = dt0.strftime(time_format)
# result_dict['var_end_time'] = dt1.strftime(time_format)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
