# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is just for adding dye to an existing file.

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
import netCDF4 as nc
import shutil
from datetime import datetime
start_time = datetime.now()

# copy rivers.nc file from riv2 to riv2dye
out_frc = Ldir['frc']
in_frc = out_frc.replace('dye','')

nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
for nc_fn in nc_list:
    in_fn = Ldir['LOogf'] + in_frc + '/' + nc_fn
    out_fn = Ldir['LOogf_f'] + nc_fn
    print('in_fn = ' + in_fn)
    print('out_fn = ' + out_fn)
    if os.path.exists(out_fn):
        os.remove(out_fn)
    shutil.copyfile(in_fn, out_fn)

# edit just one file
out_fn = Ldir['LOogf_f'] + 'ocean_clm.nc'
ds = nc.Dataset(out_fn, 'a')
vn = 'dye_01'
if vn not in ds.variables:
    vv = ds.createVariable(vn, float, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'dye'
    vv.time = "ocean_time"
    vv.units = "kg m-3"
    vv[:] = 0 * ds['salt'][:]
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

# ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)


