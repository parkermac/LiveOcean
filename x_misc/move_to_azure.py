#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

** Generic tool to copy a single file to Azure. **

Example of how to run from the linux command line:

python move_to_azure.py -fn [full path to filename]

Screen output gives the link to use for downolading the file.

With keyword arguments you can change:
-fn the input filename (typically the full path)
-fn_out the output filename (the default is to create the input filename without the path)
-c the container name (no dots or underscores allowed)
All are optional but to actually do something you need to provide the input filename

To test it try:
python move_to_azure.py
- will move the testfile to azure and return its URL to the screen
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
input_filename = Ldir['LO'] + 'extract/az_testfile.txt'
output_filename = ''
container_name = 'pm-share'

parser = argparse.ArgumentParser()
# optional command line arguments, can be input in any order, or omitted
parser.add_argument('-fn', '--input_filename', nargs='?', type=str, default=input_filename)
parser.add_argument('-fn_out', '--output_filename', nargs='?', type=str, default=output_filename)
parser.add_argument('-c', '--container_name', nargs='?', type=str, default=container_name)
args = parser.parse_args()

container_name = args.container_name
input_filename = args.input_filename
output_filename = args.output_filename

if len(output_filename) == 0:
    output_filename = input_filename.split('/')[-1]
else:
    pass # use the user input

# push results to Azure

print('\n*** Pushing selected file to Azure ***\n')
print(' - path to input file: ' + input_filename)
print(' - output filename: ' + output_filename)

# Azure commands
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
blob_service = BlockBlobService(account_name=account, account_key=key)
blob_service.create_container
blob_service.set_container_acl(container_name, public_access=PublicAccess.Container)

# write it to Azure
try:
    bname = open(input_filename, 'rb')
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

