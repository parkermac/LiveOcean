#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

** Generic tool to copy a single file to Azure. **

Example of how to run from the linux command line:

python move_to_azure.py

- will move the testfile to azure and return its URL to the screen

With keyword arguments you can change:

- the container name (no dots or underscores allowed)

- the input filename (need the full path)
 
- the output filename (the default is to use the input filename without the path)

All are optional but to actually do something you need to provide the input filename:

python move_to_azure.py -f2 [full path to filename]

Screen output gives the link to use for downolading the file.

"""

import argparse

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

# Command line arguments

# set defaults
container_name = 'pm-share'
path_to_input_file = Ldir['LO'] + 'extract/az_testfile.txt'
output_filename = ''

# optional command line arguments, can be input in any order, or omitted
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--container_name', nargs='?', type=str, default=container_name)
parser.add_argument('-f1', '--path_to_input_file', nargs='?', type=str, default=path_to_input_file)
parser.add_argument('-f2', '--output_filename', nargs='?', type=str, default=output_filename)
args = parser.parse_args()

container_name = args.container_name
path_to_input_file = args.path_to_input_file
output_filename = args.output_filename

if len(output_filename) == 0:
    output_filename = path_to_input_file.split('/')[-1]
else:
    pass # use the user input

# push results to Azure

print('\n*** Pushing selected file to Azure ***\n')
print(' - input file: ' + path_to_input_file)
print(' - output filename: ' + output_filename)

# Azure commands
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
blob_service = BlockBlobService(account_name=account, account_key=key)
blob_service.create_container(container_name)
blob_service.set_container_acl(container_name, public_access=PublicAccess.Container)

# write it to Azure
try:
    bname = open(path_to_input_file, 'rb')
    blob_service.create_blob_from_stream(container_name, output_filename, bname)
    bname.close()
    print('\n*** success ***')
    az_url = ('https://pm2.blob.core.windows.net/'
        + container_name + '/' + output_filename)
    print('\nUSE THIS URL TO ACCESS THE FILE\n')
    print(az_url)
except:
    # could be FileNotFoundError from open, or an Azure error
    print('*** FAILED TO WRITE ***')
print('\n' + 50*'*')

