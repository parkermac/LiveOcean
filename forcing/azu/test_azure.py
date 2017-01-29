"""
Azure test code.

Run from the python command line with a command like:
    
run /Users/PM5/Documents/LiveOcean/forcing/azu/test_azure.py -d 2015.04.22

"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# **************************************************

from azure.storage.blob import BlobService

#azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
#account = azu_dict['account']
#key = azu_dict['key']

f_string = 'f' + Ldir['date_string']
ff_string = f_string.replace('.','')

#containername = ff_string
#
## get a handle to your account
#blob_service = BlobService(account_name=account, account_key=key)
## When you create the container it works even if the container exists.
## I suspect it does not clobber existing files, so that could be a problem
## in some instances.
#blob_service.create_container(containername)
#blob_service.set_container_acl(containername, x_ms_blob_public_access='container')

# account name and key
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']

containername = ff_string

# get a handle to your account
blob_service = BlobService(account_name=account, account_key=key)
blob_service.create_container(containername)
blob_service.set_container_acl(containername, x_ms_blob_public_access='container')

#list all blobs of the container
blobs = blob_service.list_blobs(containername)
for blob in blobs:
    print(blob.name)
    print(blob.url)
