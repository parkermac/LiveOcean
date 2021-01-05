"""
This is the main program for running daily particle tracking jobs,
and then processing them for the interactive DRIFTER website pages.

Testing on mac:

run make_forcing_main.py -d 2019.07.04 -test True

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

# imports
from datetime import datetime, timedelta
from time import time, sleep
import subprocess
import numpy as np
import json
import netCDF4 as nc

start_time = datetime.now()
print('*** Creating drifter files for ' + Ldir['date_string'] + ' ***')
if Ldir['testing']:
    dtt = 1
    dsr = Ldir['date_string']
    if False:
        # check on days that the real run would do
        dtr00 = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
        dtr0 = dtr00 - timedelta(days = 11)
        dtr1 = dtr0 + timedelta(days = 13)
        print(dtr0.strftime('%Y.%m.%d'))
        print(dtr1.strftime('%Y.%m.%d'))
        sys.exit()
else:
    dtt = 14
    dtr00 = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    dtr0 = dtr00 - timedelta(days = 11)
    dsr = dtr0.strftime('%Y.%m.%d')

# RUN TRACKER JOBS - parallelize using subprocess
os.chdir(Ldir['LO'] + 'tracker2')
start = time()
procs = []
exp_list = ['full', 'PS']
# function defining the subprocess
def run_sub(exp, dtt):
    cmd = ['python', 'tracker.py', '-exp', exp, '-ds', dsr, '-dtt', str(dtt),
        '-clb', 'True', '-t', 'forWeb']
    proc = subprocess.Popen(cmd)
    return proc
for exp in exp_list:
        sleep(1) # cludge: needed so that calls to Lstart() don't collide while writing lo_info.csv
        proc = run_sub(exp, dtt)
        procs.append(proc)
# the proc.communicate() method will only return after a job is done
# so this loop effectively checks on all jobs sequentially, and does not
# end until they all are
for proc in procs:
    proc.communicate()
end = time()
print ('Finished tracker runs in %0.1f sec' % (end - start))
os.chdir(Ldir['LO'] + 'forcing/drifters')
out_fn_dict = {}
out_json_dict = {}
for exp in exp_list:
    out_fn_dict[exp] = (Ldir['LOo'] + 'tracks2/' + exp + '_surf_forWeb/release_' + dsr + '.nc')
    out_json_dict[exp] = (Ldir['LOo'] + 'tracks2/' + exp + '_surf_forWeb/tracks_' + exp + '.json')
result = 'success'
for exp in exp_list:
    if not os.path.isfile(out_fn_dict[exp]):
            result = 'fail'
# END OF TRACKER JOBS

# CONVERT TO JSON AND SCP TO HOMER
for exp in exp_list:
    if exp == 'full':
        skp = 3
    else:
        skp = 2
    fn = out_fn_dict[exp]
    ds = nc.Dataset(fn)
    # packed time, particle
    x = ds['lon'][:]
    y = ds['lat'][:]
    NT, NP = x.shape
    xy = []
    for pp in range(NP):
        xy.append({'x': list(x[::skp,pp]), 'y': list(y[::skp,pp])})
    json.dump(xy, open(out_json_dict[exp], 'w'))
    # send to homer
    cmd_list = ['scp',out_json_dict[exp],
        'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/tracks/'+ out_json_dict[exp].split('/')[-1]]
    ret = subprocess.call(cmd_list)
    if ret != 0:
        print('WARNING: problem moving json to homer ' + out_json_dict[exp].split('/')[-1])
        result = 'fail'
# END CONVERT AND SCP

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
result_dict['result'] = result

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
    

