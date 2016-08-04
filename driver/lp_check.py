# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:16:57 2016

@author: PM5

Check for missing low_pass files.

"""

#%% setup
import os
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
parser.add_argument("-x", "--ex_name", type=str, default='lo1')
args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

#%% look for all forecast days

dt0 = datetime(2013,1,1)
dt1 = datetime.now()

dt = dt0
while dt <= dt1:
    f_string = 'f' + datetime.strftime(dt,'%Y.%m.%d')
    if os.path.isfile(Ldir['roms'] + 'output/' + Ldir['gtagex'] +
    '/' + f_string + '/low_pass.nc'):
        pass
    else:
        print(f_string)
    dt = dt + timedelta(days=1)

