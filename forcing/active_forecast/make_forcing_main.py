#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for creating movies and images from a forecast

For testing on my mac run in ipython as

run make_forcing_main.py -g cas4 -t v2 -x lo6biom -d 2018.09.29
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
from time import time

container_name = 'active-forecast'

print(' - Creating wesite images for ' + Ldir['date_string'])
os.chdir(Ldir['LO'] + 'plotting/')

# NOTE 2020.01.24: P_tracks_MERHAB is just for Ryan, it does not go to the website
#P_list = ['P_3day', 'P_tracks_MERHAB', 'P_merhab2', 'P_tracks_ps', 'P_willapa_omega']
P_list_1 = ['P_tracks_MERHAB', 'P_merhab2', 'P_willapa_omega']
P_list_2 = ['P_basic', 'P_Chl_DO', 'P_basic_salish', 'P_Chl_DO_salish', 'P_tracks_barber']


P_list = P_list_1 + P_list_2
#P_list = P_list_2
#P_list = ['P_basic','P_Chl_DO']

result = 'success'
for P_name in P_list:
    
    tt0 = time()
    
    print('\n  -- working on %s' % (P_name))
    
    if P_name in ['P_tracks_MERHAB', 'P_merhab2', 'P_tracks_barber']:
        lt = 'merhab'
    else:
        lt = 'forecast'
        
    if P_name == 'P_3day':
        do_mov = 'False'
    else:
        do_mov = 'True'
        
    cmd = ['python','pan_plot.py',
                  '-g', Ldir['gridname'],'-t', Ldir['tag'],'-x', Ldir['ex_name'],
                  '-lt', lt,'-0', Ldir['date_string'],'-1', Ldir['date_string'],
                  '-pt', P_name,'-mov', do_mov]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if False:
        print('\n-main: screen output from subprocess-')
        print(proc.stdout.decode())
        print('\n-main: errors from subprocess-')
        # for some reason the ffmpeg output ends up in stderr
        print(proc.stderr.decode())
    
    if P_name == 'P_3day':
        fn = Ldir['LOo'] + 'plots/' + lt + '_' + P_name + '_' + Ldir['gtagex'] + '/plot_0000.png'
        out_fn = P_name + '.png'
    else:
        fn = Ldir['LOo'] + 'plots/' + lt + '_' + P_name + '_' + Ldir['gtagex'] + '/movie.mp4'
        out_fn = P_name + '.mp4'
    
    #result = write_to_azure(fn, blob_service, containername, out_fn)
    
    input_filename = fn #in_dir + output_filename
    output_filename = out_fn
    
    az_dict = Lfun.copy_to_azure(input_filename, output_filename, container_name, Ldir)
    if az_dict['result'] == 'fail':
        result = 'fail'
        print('Failed to copy to Azure:')
        print(output_filename)
    
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

