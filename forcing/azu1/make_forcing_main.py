#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the main program for pushing files to AZURE.

For testing on my mac run in ipython as

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

print(' - Pushing selected files to Azure for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']
ff_string = f_string.replace('.','') # azure does not like dots in container names
container_name = ff_string

# input directory
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'

# files to copy to azure
file_list = ['ocean_surface.nc', 'low_passed_UBC.nc']

result = 'success'
for output_filename in file_list:
    input_filename = in_dir + output_filename
    az_dict = Lfun.copy_to_azure(input_filename, output_filename, container_name, Ldir)
    if az_dict['result'] == 'fail':
        result = 'fail'
        print('Failed to copy to Azure:')
        print(output_filename)
        
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

# ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)

