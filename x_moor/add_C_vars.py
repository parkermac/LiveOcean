# -*- coding: utf-8 -*-
"""
Add carbon variables to a single mooring file.
"""


import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

# choose file to process
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
infile = Lfun.choose_item(indir, tag='.nc')

#%% calculate pH and Aragonite saturation state
import subprocess
func = ("run_co2sys(\'" +
    indir + "\',\'" +
    infile + "\',\'" +
    'input3' + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
print(out.decode())
print(err.decode())
