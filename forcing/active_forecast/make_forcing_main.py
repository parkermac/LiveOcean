#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for creating movies and images from a forecast

To test on boiler from ipython:

run make_forcing_main.py -d 2020.08.10
"""

import os, sys
sys.path.append(os.path.abspath('../'))
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

testing = True

from datetime import datetime, timedelta
start_time = datetime.now()
import subprocess
from time import time

print(' - Creating wesite images for ' + Ldir['date_string'])
os.chdir(Ldir['LO'] + 'plotting/')

# NOTE 2020.01.24: P_tracks_MERHAB is just for Ryan, it does not go to the website
P_list_1 = ['P_tracks_MERHAB', 'P_merhab2', 'P_willapa_omega']
P_list_2 = ['P_basic', 'P_Chl_DO', 'P_basic_salish', 'P_Chl_DO_salish', 'P_tracks_barber']
P_list = P_list_1 + P_list_2

if testing:
    P_list = ['P_DO']

result = 'success'
for P_name in P_list:
    tt0 = time()
    print('\n  -- working on %s' % (P_name))
    
    dstr0 = Ldir['date_string']
    dstr1 = Ldir['date_string']
    
    if P_name in ['P_tracks_MERHAB', 'P_merhab2', 'P_tracks_barber']:
        lt = 'merhab'
    elif P_name == 'P_DO':
        lt = 'daily_plus'
        dt1 = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
        dt0 = dt1 - timedelta(days=60)
        dstr0 = dt0.strftime(format='%Y.%m.%d')
    else:
        lt = 'forecast'
        
    # if testing:
    #     lt = 'allhours'
    
    # run the plotting code
    cmd = ['python','pan_plot.py',
                  '-g', Ldir['gridname'],'-t', Ldir['tag'],'-x', Ldir['ex_name'],
                  '-lt', lt,'-0', dstr0,'-1', dstr1,
                  '-pt', P_name,'-mov', 'True']
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if False:
        print('\n-main: screen output from subprocess-')
        print(proc.stdout.decode())
        print('\n-main: errors from subprocess-')
        # for some reason the ffmpeg output ends up in stderr
        print(proc.stderr.decode())
    
    input_filename = Ldir['LOo'] + 'plots/' + lt + '_' + P_name + '_' + Ldir['gtagex'] + '/movie.mp4'
    output_filename = P_name + '.mp4'
    
    if P_name == 'P_tracks_MERHAB':
        # send file to azure
        print('\nCopying '+output_filename+' to azure')
        container_name = 'active-forecast'
        az_dict = Lfun.copy_to_azure(input_filename, output_filename, container_name, Ldir)
        if az_dict['result'] == 'fail':
            result = 'fail'
            print('Failed to copy to Azure:')
            print(output_filename)
    else:
        # send file to homer (only works from boiler)
        print('\nCopying '+output_filename+' to homer')
        cmd_list = ['scp',input_filename,
            'pmacc@homer.u.washington.edu:/hw00/d47/pmacc/LO/Figs_active_forecast/'+output_filename]
        ret = subprocess.call(cmd_list)
        print('Return code = ' + str(ret) + ' (0=success)')
    
    print('  -- took %0.1f seconds' % (tt0-time()))
    
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

