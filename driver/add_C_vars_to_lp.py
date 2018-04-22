#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 15:42:30 2018

@author: pm7

Add PH and ARAG to low_passed files
"""

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str,
                    default='cascadia1')
parser.add_argument('-t', '--tag', nargs='?', type=str,
                    default='base')
parser.add_argument('-x', '--ex_name', nargs='?', type=str,
                    default='lobio5')
parser.add_argument('-d', '--date_string', nargs='?', type=str,
                    default='2017.01.15')

#from importlib import reload
#reload(Lfun)

Ldir = Lfun.Lstart()

args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + args.date_string + '/'

fn = indir + 'low_passed.nc'

#%% calculate pH and Aragonite saturation state
import subprocess
func = ("run_co2sys(\'" +
    indir + "\',\'" +
    'low_passed.nc' + "\',\'" +
    'input3' + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
print(out.decode())
print(err.decode())

#zfun.ncd(fn)