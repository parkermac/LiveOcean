#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 07:52:18 2017

@author: PM5

Plot raw and filtered fields from the hycom1 archive.
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd

import matplotlib.pyplot as plt

# specify the output directory
out_dir = Ldir['data'] + 'hycom1/'

k0 = 5
df = pd.read_pickle(out_dir + 'extraction_raw_' + str(k0) + '.p')
df_fh = pd.read_pickle(out_dir + 'extraction_fh_' + str(k0) + '.p')

t0 = df_fh.index[0]
t1 = df_fh.index[-1]

#%% combine
for vn in df_fh.keys():
    df[vn + '_fh'] = df_fh[vn]
    
#%% plotting

plt.close('all')

fig, axes = plt.subplots(5, 1, sharex=True)

nr = 0
for vn in df_fh.keys():
    ax = axes[nr]
    df.ix[:,[vn, vn+'_fh']].plot(ax=ax, grid=True, xlim=(t0, t1))
    nr+=1



