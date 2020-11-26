"""
This is the main program for making the CRITFC output file,
using code from Charles Seaton.

For testing on my mac run in ipython as
run make_forcing_main.py -d 2019.07.04 -test True

"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

import subprocess
from datetime import datetime
import zrfun

start_time = datetime.now()

out_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + Ldir['date_string'] + '/'

# critfc code inputs
rundate = Ldir['date_string'].replace('.','-')
depthfile = out_dir + 'ocean_his_0001.nc'
hgrid = Ldir['data'] + 'critfc/hgrid.ll'
vgrid = './vgrid.in'
basedir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
outdir = out_dir

cmd = ['python', 'gen_cmop_nudge.py', hgrid, vgrid, depthfile,
    basedir, outdir, rundate,'-test',str(Ldir['testing'])]
proc = subprocess.Popen(cmd)
proc.communicate()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if proc.returncode == 0:
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)