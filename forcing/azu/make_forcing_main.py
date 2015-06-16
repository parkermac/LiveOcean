"""
This is the main program for making the AZU forcing file.
"""
# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("frc", type=str, help="atm, ocn, riv, or tide")
parser.add_argument("run_type", type=str, help="forecast or backfill")
parser.add_argument("date_string", type=str, help="e.g. 2014.02.14")
args = parser.parse_args()

# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
Ldir['LOogf_f'] = (Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/' + args.frc + '/')
    
# screen output
from datetime import datetime
print('MAIN: frc = ' + args.frc + ', run_type = ' + args.run_type
    + ', date_string = ' + args.date_string)
print('MAIN start time = ' + str(datetime.now()))

# ****************** CASE-SPECIFIC CODE *****************

# Azure commands
import azure.storage as az    
print('\nPushing forecast files to Azure for ' + args.date_string + '\n')
# azure does not like dots in container names!
f_string = 'f' + args.date_string
ff_string = f_string.replace('.','')

# account name and key
#azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_rob_2015.04.24.csv')
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']

containername = ff_string

# get a handle to your account
blob_service = az.BlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, x_ms_blob_public_access='container')

if False: # testing
    nend = 3
else:
    if run_type == 'backfill':
        nend = 26
    elif run_type == 'forecast':
        nend = 74 # new three-day forecast 5/25/2015
    
result_dict = dict()
try:        
    for ii in range(2,nend):
        ncpad = '0000' + str(ii)
        ncpad = ncpad[-4:]
        hisname = 'ocean_his_' + ncpad + '.nc'
        dirname = Ldir['roms'] + 'output/' + Ldir['gtag'] + '/' + f_string + '/'
        fname = dirname + hisname   
        bname = open(fname, 'r')
        blob_service.put_block_blob_from_file(containername, hisname, bname)
        print 'done putting ' + hisname
        bname.close()
    result_dict['result'] = 'success'
except:
    result_dict['result'] = 'fail'
    
# write results to an output file for the driver
csv_name_out = Ldir['LOogf_f'] + 'Info/' + 'process_status.csv' 
Lfun.dict_to_csv(result_dict, csv_name_out)

# ******************************************************* 
# no worker code in this case
print('MAIN end time = ' + str(datetime.now()))

