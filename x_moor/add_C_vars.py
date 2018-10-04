# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:28:25 2016

@author: PM5

Add carbon variables to a mooring file.
"""

# setup
#import netCDF4 as nc
#import matplotlib.pyplot as plt
#import numpy as np
#import matplotlib.dates as mdates
#from datetime import datetime

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
#from importlib import reload
#reload(Lfun)

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
# choose the mooring extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
infile = Lfun.choose_item(indir, tag='.nc')
fn = indir + infile

#%% calculate pH and Aragonite saturation state
import subprocess
func = ("run_co2sys(\'" +
    indir + "\',\'" +
    infile + "\',\'" +
    'input3' + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
print(out.decode())
print(err.decode())

import netCDF4 as nc
ds = nc.Dataset(fn)
for vn in ds.variables:
    print(vn)
ds.close()