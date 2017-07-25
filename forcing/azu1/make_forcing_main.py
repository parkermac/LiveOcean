#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:34:50 2017

@author: PM5

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

# Azure commands
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess
ff_string = f_string.replace('.','') # azure does not like dots in container names
# account name and key
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
containername = ff_string
# get a handle to the account
blob_service = BlockBlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, public_access=PublicAccess.Container)

# input directory
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
# output files
out_list = ['ocean_surface.nc', 'low_passed.nc', 'low_passed_UBC.nc']

in_dir2 = Ldir['LOo'] + 'plots/merhab_P_tracks_MERHAB_' + Ldir['gtagex'] + '/'
out_list2 = ['movie.mp4']

def write_to_azure(out_fn, blob_service, containername, outname):
    # write it to Azure
    try:
        bname = open(out_fn, 'rb')
        blob_service.create_blob_from_stream(containername, out_name, bname)
        print('done putting ' + out_name)
        bname.close()
    except:
        # could be FileNotFoundError from open, or an Azure error
        print(' - Unable to write ' + out_name + ' to Azure')
        result = 'fail'
    return result

result = 'success'

for out_name in out_list:
    out_fn = in_dir + out_name
    result = write_to_azure(out_fn, blob_service, containername, out_name)
    
for out_name in out_list2:
    out_fn = in_dir2 + out_name
    result = write_to_azure(out_fn, blob_service, containername, out_name)
        
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

