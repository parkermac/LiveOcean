#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 13:58:44 2018

@author: pm7

Code to copy all ubc files into one place, renaming them by
forecast day.
"""


#%% setup
import os
import shutil
import sys
import argparse
from datetime import datetime, timedelta

alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

# get command line arguments, if any
parser = argparse.ArgumentParser()
# optional arguments
parser.add_argument("-g", "--gridname", type=str, default='cascadia1')
parser.add_argument("-t", "--tag", type=str, default='base')
parser.add_argument("-x", "--ex_name", type=str, default='lobio5')
args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# create output directory
out_dir = Ldir['LOo'] + 'UBC_all_' + Ldir['gtagex'] + '/'
Lfun.make_dir(out_dir, clean=True)

# look for all forecast days

if True:
    dt0 = datetime(2013,1,2)
    dt1 = datetime.now()
else:
    dt0 = datetime(2013,1,2)
    dt1 = datetime(2013,1,5)

dt = dt0
while dt <= dt1:
    f_string = 'f' + datetime.strftime(dt,'%Y.%m.%d')
    fdir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
    in_fn = fdir + 'low_passed_UBC.nc'
    out_fn = out_dir + 'low_passed_UBC_' + f_string + '.nc'
    if os.path.isfile(in_fn):
        shutil.copyfile(in_fn, out_fn)
        print('copied file for ' + f_string)
    else:
        print(' - MISSING file for ' + f_string)
        
    dt = dt + timedelta(days=1)

