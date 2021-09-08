#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for creating movies and images from a forecast

To test on boiler from ipython:
run make_forcing_main.py -d [YYYY.MM.DD for today] -test True

To run on boiler for today's forecast, from command line:
python make_forcing_main.py -d [YYYY.MM.DD for today] > log &

To test on mac from ipython:
cd /Users/pm8/Documents/LiveOcean/forcing/active_forecast2
run make_forcing_main.py -d 2019.07.04 -test True

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

ds00 = Ldir['date_string']
dt00 = datetime.strptime(ds00,'%Y.%m.%d')
    
procs = []

if False:#Ldir['testing'] == True:
    moviename_list = ['P1_nshelf_oxygen_bot']
else:
    moviename_list = ['P1_full_salt_top', 'P1_full_oxygen_bot', 'P1_nshelf_oxygen_bot',
                        'P1_PS_temp_top', 'P1_PS_speed_top',
                        'P1_willapa_ARAG_top', 'P1_willapa_ARAG_bot',
                        'Phab_full_salt_top']
                    
os.chdir(Ldir['LO'] + 'plot5/')

if Ldir['testing'] == False:
    tt0 = time()
    result = 'success'
    for moviename in moviename_list:
        (pt,dom,vn,BOT) = moviename.split('_')
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
        
        if moviename == 'Phab_full_salt_top':
            ttag = 'hab'
        else:
            ttag = 'base'
        
        if moviename == 'P1_nshelf_oxygen_bot':
            # start this one a couple day s earlier
            dt0 = dt00 - timedelta(days=2) # start earlier
            dt1 = dt00 + timedelta(days=Ldir['forecast_days']-1)
        else:
            dt0 = dt00
            dt1 = dt0 + timedelta(days=Ldir['forecast_days']-1)
        ds0 = dt0.strftime('%Y.%m.%d')
        ds1 = dt1.strftime('%Y.%m.%d')
        
        cmd = ['python', 'p5.py', '-ds0', ds0, '-ds1', ds1, '-lt', 'hourly', '-mov', 'True',
            '-pt', pt,
            '-dom', dom, '-vn', vn, '-tracks', tracks, '-emask', emask, '-ttag', ttag,
            '-avl', avl, '-bot', bot]
        
        print('\n' + moviename)
        sys.stdout.flush()
        sleep(4)
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

    if True:#Ldir['testing'] == False:
        sleep(4)
        # send file to homer (only works from boiler)
        print('\nCopying '+output_filename+' to homer')
        sys.stdout.flush()
        try:
            cmd_list = ['scp',input_filename,
                'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast/'+output_filename]
                
            proc = subprocess.Popen(cmd_list,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # note "error" appears to be where the ffmpeg screen output goes
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                print(' sdtout '.center(60,'-'))
                print(stdout.decode())
            if len(stderr) > 0:
                print(' stderr '.center(60,'-'))
                print(stderr.decode())
            sys.stdout.flush()
            # ret = subprocess.call(cmd_list)
            # if ret != 0:
            #     print('WARNING: problem with movie ' + moviename)
            #     result = 'fail'
            # print('Return code = ' + str(ret) + ' (0=success)')
        except Exception as e:
            print(e)
            pass
    
    # and save a local copy
    print(' -- saving local copy of movie')
    try:
        shutil.copyfile(input_filename, alt_outdir + output_filename)
    except Exception as e:
        print(e)
        pass

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

