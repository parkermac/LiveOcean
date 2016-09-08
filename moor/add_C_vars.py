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
#from importlib import reload
#reload(Lfun)

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to add C vars to **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if '.nc' in m:
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
moor_file = m_dict[my_npt]
fn = indir + moor_file

#%% calculate pH and Aragonite saturation state
import subprocess
func = ("run_co2sys(\'" +
    indir + "\',\'" +
    moor_file + "\',\'" +
    'input3' + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
#print(out.decode())

#zfun.ncd(fn)