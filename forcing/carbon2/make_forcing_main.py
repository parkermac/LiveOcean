"""
This is the main program for making the CARBON variable additions
to the ROMS history files.  Designed to use the new co2sys programs.

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

pth = os.path.abspath('../co2sys')
if pth not in sys.path:
    sys.path.append(pth)
import carbon_fun as cfun
from importlib import reload
reload(cfun)

from datetime import datetime
start_time = datetime.now()

import subprocess

Ldir['indir'] = (Ldir['roms'] + 'output/' + Ldir['gtagex']
            + '/f' + Ldir['date_string'] + '/')

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

# make a list of tuples of (start,end) history file numbers
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

# --------- do the calc ---------------------------------

from multiprocessing.pool import ThreadPool
import subprocess

def work(sample):
    my_tool_subprocess = subprocess.Popen('mytool {}'.format(sample),shell=True, stdout=subprocess.PIPE)
    line = True
    while line:
        myline = my_tool_subprocess.stdout.readline()
        #here I parse stdout..


num = None  # set to the number of workers you want (it defaults to the cpu count of your machine)
tp = ThreadPool(num)
for sample in all_samples:
    tp.apply_async(work, (sample,))

tp.close()
tp.join()

# (1) with a single fn in hand, and some counter "NN" that can
# be used to safely put the temp output in, without being overwritten
# by another subprocess
vv_dict = cfun.get_full(fn)
tempdir0 = Ldir['LOo'] + 'co2sys/'
Lfun.make_dir(tempdir0)
tempdir = tempdir0 + 'temp_' + str(NN) + '/'
Lfun.make_dir(tempdir, clean=True) # NOTE: make_carbon() does this, no nead to do twice
cfun.make_carbon(vv_dict, tempdir, print_info=True)
PH, OM = cfun.get_carbon(tempdir)

# (2) save the output to fn
ds = nc.Dataset(fn, 'a')
if


# -------------------------------------------------------

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

if Ldir['run_type'] == 'low_passed':
    ds0 = nc.Dataset(Ldir['indir'] + 'low_passed.nc')
    if ('PH' in ds0.variables):
        result_dict['result'] = 'success'   
    else:
        result_dict['result'] = 'fail' 
    ds0.close()
else:    
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
