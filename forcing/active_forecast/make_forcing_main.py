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

# Azure commands
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess
# account name and key
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
containername = 'active-forecast'
# get a handle to the account
blob_service = BlockBlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, public_access=PublicAccess.Container)

def write_to_azure(out_fn, blob_service, containername, out_name):
    # write it to Azure
    try:
        bname = open(out_fn, 'rb')
        blob_service.create_blob_from_stream(containername, out_name, bname)
        print('done putting ' + out_name)
        bname.close()
        result = 'success'
    except:
        # could be FileNotFoundError from open, or an Azure error
        print(' - Unable to write ' + out_name + ' to Azure')
        result = 'fail'
    return result

print(' - Creating wesite images for ' + Ldir['date_string'])
os.chdir(Ldir['LO'] + 'plotting/')

P_list = ['P_3day', 'P_tracks_MERHAB', 'P_merhab2', 'P_tracks_ps', 'P_willapa_omega', 'P_rho']
#P_list = ['P_willapa_omega']
#P_list = ['P_tracks_ps']
#P_list = ['P_3day']
#P_list = ['P_rho']
for P_name in P_list:
    
    if P_name in ['P_tracks_MERHAB', 'P_merhab2', 'P_tracks_ps']:
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
        
    result = write_to_azure(fn, blob_service, containername, out_fn)
    
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

