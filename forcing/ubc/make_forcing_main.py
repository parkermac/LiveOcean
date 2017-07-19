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
from datetime import datetime, timedelta
start_time = datetime.now()

# make input name (full path)
in_fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
    '/f' + Ldir['date_string'] + '/low_passed.nc')
# and also the output file name
out_fn = in_fn.replace('.nc', '_UBC.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# do the extraction (not that you have to feed it a list, even for one item)
import UBC_subdomain
UBC_subdomain.get_UBC_subdomain([in_fn,])

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