#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for creating the MERHAB TRACKS plots and movie.

For testing on my mac run in ipython as

cd ~/Documents/LiveOcean/forcing/tracks_m

run make_forcing_main.py -d 2017.05.18
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

print(' - Creating MERAB movie for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']

os.chdir(Ldir['LO'] + 'plot1/')

cmd = ['python','pan_plot.py',
              '-g', Ldir['gridname'],
              '-t', Ldir['tag'],
              '-x', Ldir['ex_name'],
              '-lt', 'merhab',
              '-0', Ldir['date_string'],
              '-pt', 'P_tracks_MERHAB',
              '-avl', 'False',
              '-mov', 'True']
              
proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print('\n-main: screen output from subprocess-')
print(proc.stdout.decode())
print('\n-main: errors from subprocess-')
# for some reason the ffmpeg output ends up in stderr
print(proc.stderr.decode())

# move the movie and the last frame to the forecast folder
in_dir = Ldir['LOo'] + 'plots/merhab_P_tracks_MERHAB_' + Ldir['gtagex'] + '/'
out_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
fn_list = ['movie.mp4', 'plot_0070.png']
for fn in fn_list:
    cmd2 = ['cp',in_dir+fn,out_dir+fn]
    proc2 = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if proc2.returncode == 0:
    result = 'success'
else:
    result = 'fail'
    
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

