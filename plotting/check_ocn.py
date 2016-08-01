# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 07:14:06 2016

@author: PM5

Code to look at the ocn files on a given day.

"""

#%% setup

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import zfun; reload(zfun)

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

#%% where to look

dir0 = '/Users/PM5/Documents/LiveOcean_output/'
case = 1
if case == 1:
    gtag = 'cascadia2_frc2'
    f_string = 'f2013.01.01'
elif case == 2:
    gtag = 'cascadia1_base'
    f_string = 'f2015.09.19'

#%% look at the output files

in_fn = (dir0 + gtag + '/' + f_string + '/ocn/ocean_clm.nc')

ds = nc.Dataset(in_fn)

#zfun.ncd(ds)

#%% plotting

plt.close()

vn_list =['salt', 'temp', 'u', 'v']

NP = len(vn_list)
NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
NC = np.ceil(np.sqrt(NP)).astype(int)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,10), squeeze=False)
cc = 0
for vn in vn_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    # **********************
    v = ds[vn][0, -1, :, :]
    print(vn + ' ' + str(np.isnan(v).sum()) + ' nans')
    cs = ax.pcolormesh(v)
    ax.axis('tight')
    fig.colorbar(cs, ax=ax)
    ax.set_title(vn)
    # **********************
    cc += 1
plt.show()



plt.show()

