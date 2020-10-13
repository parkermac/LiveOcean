"""
This is the main program for making a SURFACE subset of the daily output.

It creates a single NetCDF file containing only the surface fields
from all the history files in a given day.

Testing on mac:

run make_forcing_main.py -d 2019.07.04

2019.07.18 Shortened the list of variables to only those served by EDS Viewer,
and used every other hour.

2020.07.14 Added a second output file for SCOOT, including layers at selected z-levels.
Also cleaned up code and created module surf_fun.py.
"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

# imports
import surf_fun
from datetime import datetime
import netCDF4 as nc
import zrfun
import numpy as np
sys.path.append(os.path.abspath('../../plotting/'))
import pfun
from time import time

start_time = datetime.now()

print(' - Creating surface file(s) for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']

# Create out_dir
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
out_dir = in_dir # output goes to same place as input

# ======== Create an output file for EDS ===============================
fn_list_raw = os.listdir(in_dir)
fn_list = []
for item in fn_list_raw:
    if 'ocean_his' in item and '.nc' in item:
        fn_list.append(in_dir + item)
fn_list.sort()
# shorten the list to be every 2 hours
fn_list = fn_list[::2]

# Initialize the multi-file input dataset
in_ds = nc.MFDataset(fn_list)

# Initialize
out_fn = out_dir + 'ocean_surface.nc'
print(' - Writing to: ' + out_fn)
out_ds = surf_fun.create_ds(out_fn)
surf_fun.add_dims(in_ds, out_ds)
# Add standard 3D-T, rho-grid, variables
vn_list2t = ['Uwind', 'Vwind', 'ocean_time']
vn_list3t = ['salt', 'temp']
slev = -1
surf_fun.add_fields(in_ds, out_ds, vn_list2t, vn_list3t, slev=slev)
# Add custom variables
# - surface velocity on the rho grid
u0 = in_ds['u'][:, slev, :, :].squeeze()
v0 = in_ds['v'][:, slev, :, :].squeeze()
u = np.nan * in_ds['salt'][:, slev, :, :].squeeze()
v = u.copy()
u[:, :, 1:-1] = (u0[:, :, 1:] + u0[:, :, :-1])/2
v[:, 1:-1, :] = (v0[:, 1:, :] + v0[:, :-1, :])/2
#
vv = out_ds.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'eastward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = u
#
vv = out_ds.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'northward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = v
# Close output Dataset
out_ds.close()
# Close multi-file input Dataset
in_ds.close()
# ======================================================================

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

