#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for creating movies and images from a forecast

To test on boiler from ipython:

run make_forcing_main.py -d [YYYY.MM.DD for today]
"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime, timedelta
start_time = datetime.now()
import subprocess
from time import time, sleep
import shutil

print(' - Creating wesite images for ' + Ldir['date_string'])

ds0 = Ldir['date_string']
# for the purposes of this driver we always assume we are plotting over a forecast
dt0 = datetime.strptime(ds0,'%Y.%m.%d')
dt1 = dt0 + timedelta(days=Ldir['forecast_days']-1)
ds1 = dt1.strftime('%Y.%m.%d')
    
procs = []
moviename_list = ['full_salt_top', 'full_oxygen_bot',
                    'PS_temp_top', 'PS_speed_top',
                    'willapa_ARAG_top', 'willapa_ARAG_bot']
                    
os.chdir(Ldir['LO'] + 'plot5/')

tt0 = time()
result = 'success'
for moviename in moviename_list:
    (dom,vn,BOT) = moviename.split('_')
    tracks = 'False'
    emask = 'False'
    avl = 'True'
    bot = 'False'
    if BOT == 'bot':
        bot = 'True'
    if vn in ['oxygen', 'ARAG', 'speed']:
        avl = 'False'
    if (vn in ['oxygen', 'ARAG']) and (dom == 'full'):
        emask = 'True'
    if vn == 'salt':
        tracks = 'True'
        
    cmd = ['python', 'p5.py', '-ds0', ds0, '-ds1', ds1, '-lt', 'hourly', '-mov', 'True',
        '-dom', dom, '-vn', vn, '-tracks', tracks, '-emask', emask,
        '-avl', avl, '-bot', bot]
        
    print('\n' + moviename)
    sys.stdout.flush()
    sleep(1)
    proc = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    procs.append(proc)

for proc in procs:
    # note "error" appears to be where the ffmpeg screen output goes
    output, errors = proc.communicate()
    
print('time to run all jobs = %0.1f sec' % (time() - tt0))
sys.stdout.flush()

alt_outdir = Ldir['LOo'] + 'Figs_active_forecast/'
Lfun.make_dir(alt_outdir, clean=True)

result = 'success'
for moviename in moviename_list:
    input_filename = Ldir['LOo'] + 'p5/' + Ldir['gtagex'] + '/' + moviename + '/' + moviename + '.mp4'
    output_filename = moviename + '.mp4'

    # send file to homer (only works from boiler)
    # print('\nCopying '+output_filename+' to homer')
    sys.stdout.flush()
    cmd_list = ['scp',input_filename,
        'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast/'+output_filename]
    ret = subprocess.call(cmd_list)
    if ret != 0:
        print('WARNING: problem with movie ' + moviename)
        result = 'fail'
    # print('Return code = ' + str(ret) + ' (0=success)')
    # print('  -- took %0.1f seconds' % (tt0-time()))
    
    # and save a local copy
    shutil.copyfile(input_filename, alt_outdir + output_filename)

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

