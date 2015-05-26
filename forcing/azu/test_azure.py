"""
Azure test code.

Run from the python command line with a command like:
    
%run /Users/PM5/Documents/LiveOcean/forcing/azu/test_azure.py 2015.04.22

"""

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("date_string", type=str, help="e.g. 2014.02.14")
args = parser.parse_args()

# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

import azure.storage as az

#azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_rob_2015.04.24.csv')
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']

f_string = 'f' + args.date_string
ff_string = f_string.replace('.','')

containername = ff_string

# get a handle to your account
blob_service = az.BlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, x_ms_blob_public_access='container')

#list all blobs of the container
blobs = blob_service.list_blobs(containername)
for blob in blobs:
    print('')
    print(blob.name)
    print(blob.url)
   
