# -*- coding: utf-8 -*-
"""
This is the main program for making the OCN forcing file.

It is designed to work with the new hycom2 archive of extracted files,
as well as regular forecasts.

Performance: on my mac a backfill day for cas6 takes 1 minute.

2017.12.10 I added a planB flag to the forecast case
which copies the clm file from the previous day and
makes the final time a day later.

2018.05.19 I added the add_CTD flag (and the Ofun_CTD module) to add CTD data
on a specified day only.

2019.04.24 I Removed the fix_NSoG code entirely because that experiment did
not produce the desired results.

2019.05.09 Changed the day for which CTD ICs are added to 2016.12.15, but in
Ofun_CTD.get_casts() it is hardwired to look for casts on or after January 2017.

2019.05.20 Added Ofun.get_interpolated_alt() which sped up the program by a factor of 10.

2020.08.13 Added a method to get the HYCOM files using the nco operator "ncks",
which is new essentially Plan A.  I deprecated the "FMRC_best" method to Plan B
because it was failing about tome time out of four.  Then moved the persistence
backup to Plan C.

2020.09.27 Moved the ncks method to Ofun, and reworked the error-handling logic.
It was working well but the fmrc method would not generate enough of an
exception to move to Plan C.

*******************************

To run from the command line in LiveOcean/driver/:
    
./driver_forcing2.sh -g cas6 -t v1 -f ocn4 -r backfill -0 20170101 -1 20170101

To test in python on mac:

# standard backfill
run make_forcing_main.py -g cas6 -t v3 -r backfill -d 2017.04.20

# backfill with Salish and coastal estuary IC's from CTD and other info
run make_forcing_main.py -g cas6 -t v1 -r backfill -d 2016.12.15
- the switch to do this is hardwired to a day: 2016.12.15

# today's forecast
run make_forcing_main.py -g cas6 -t v3 -r forecast

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


