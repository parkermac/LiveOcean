"""
This is the main program for making the CARBON variable additions
to the ROMS history files.

To run from the command line in LiveOcean/driver/:
    
./driver_forcing1.sh -g cascadia1 -t base -f carbon -r backfill -0 20170518 -1 20170518

To test in python on mac:

cd /Users/PM5/Documents/LiveOcean/forcing/carbon

run make_forcing_main.py -r backfill -d 2017.05.18

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
from datetime import datetime
start_time = datetime.now()

import subprocess

Ldir['indir'] = (Ldir['roms'] + 'output/' + Ldir['gtagex']
            + '/f' + Ldir['date_string'] + '/')

def make_func(Ldir):
    # the MATLAB command
    func = ("make_forcing_worker(\'" +
        Ldir['gridname'] + "\',\'" +
        Ldir['tag'] + "\',\'" +
        Ldir['date_string'] + "\',\'" +
        Ldir['run_type'] + "\',\'" +
        Ldir['indir'] + "\',\'" +
        Ldir['h0'] + "\',\'" +
        Ldir['h1'] + "\',\'" +
        Ldir['LOogf_f'] + "\')")
    return func

# make a list of all history files in the directory
his_list_raw = os.listdir(Ldir['indir'])
his_list = [hh for hh in his_list_raw if 'ocean_his' in hh]
his_list.sort()

# and make a list of their numbers
h_list = [int(his[-7:-3]) for his in his_list]

# set number of history files to send to a single MATLAB job
nf = 15

testing = False
if testing:
    nh = 5
    his_list = his_list[:nh]
    h_list = h_list[:nh]
    nf = 2

#%% make a list of tuples of (start,end) history file numbers
htup_list = []
ii = 0
h1 = -1 # this should always be less than h_list[-1]
while h1 < h_list[-1]:   
    if ii==0:
        h0 = h_list[0]
    else:
        h0 = h1 + 1        
    try:
        h1 = h_list[(ii+1)*nf - 1]
    except IndexError:
        h1 = h_list[-1]       
    htup_list.append((h0,h1))    
    ii += 1

# then reverse the order so that we do the first batch last
# because this will always have the longest list of items
# and so is best suited for the "waiting" branch below.
htup_list.reverse()

#%% run the subprocess
print('\n-main: start of jumbled output from first simultaneous jobs-')
for htup in htup_list:
    h0 = htup[0]
    h1 = htup[1]
    print('--output from main--')
    print('sending h0=%d : h1=%d to subprocess' % (h0, h1))
    sys.stdout.flush()
    Ldir['h0'] = str(h0)
    Ldir['h1'] = str(h1)
    func = make_func(Ldir)
    run_cmd = [Ldir['which_matlab'], "-nodisplay", "-r", func, "&"]
    if htup == htup_list[-1]:
        print(' -waiting for subprocess to complete')
        sys.stdout.flush()
        # for some reason, capturing the stdout causes the calling program
        # to wait until the subprocess finishes.
        proc = subprocess.run(run_cmd, stdout=subprocess.PIPE)
        print('\n-main: start of clean output from last simultaneous job-')
        print(proc.stdout.decode())
        sys.stdout.flush()
    else:
        subprocess.run(run_cmd)
        # whereas here the stdout goes to screen_output.txt in jumbled form
        # because several subprocesses are writing to it at once.

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)

import netCDF4 as nc
ds0 = nc.Dataset(Ldir['indir'] + his_list[0])
ds1 = nc.Dataset(Ldir['indir'] + his_list[-1])
if ('PH' in ds0.variables) and ('PH' in ds1.variables):
    result_dict['result'] = 'success'   
else:
    result_dict['result'] = 'fail' 
ds0.close()
ds1.close()

#%% ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)
