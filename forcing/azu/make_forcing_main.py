"""
This is the main program for making the AZU forcing file.
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

# Azure commands
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess

print('\nPushing forecast files to Azure for ' + Ldir['date_string'] + '\n')
# azure does not like dots in container names!
f_string = 'f' + Ldir['date_string']
ff_string = f_string.replace('.','')

# account name and key
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']

containername = ff_string

# get a handle to your account
blob_service = BlockBlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, public_access=PublicAccess.Container)

if False: # testing
    nend = 3
else:
    if Ldir['run_type'] == 'backfill':
        nend = 26
    elif Ldir['run_type'] == 'forecast':
        nend = 74 # new three-day forecast 5/25/2015

try:
    for ii in range(2,nend):
        ncpad = '0000' + str(ii)
        ncpad = ncpad[-4:]
        hisname = 'ocean_his_' + ncpad + '.nc'
        dirname = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
        fname = dirname + hisname
        bname = open(fname, 'rb')
        blob_service.put_block_blob_from_file(containername, hisname, bname)
        print('done putting ' + hisname)
        bname.close()
    result = 'success'
except:
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

