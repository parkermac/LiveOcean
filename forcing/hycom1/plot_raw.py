#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 09:08:41 2017

@author: PM5

Plot raw fields from the hycom1 archive.
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

import zfun

import pandas as pd

import matplotlib.pyplot as plt

# specify the output directory
out_dir = Ldir['data'] + 'hycom1/'

k0 = 5
df = pd.read_pickle(out_dir + 'extraction_raw_' + str(k0) + '.p')

#%% filtering

df1 = df.dropna()
t = df1.index.values
u = df1['u3d'].values

uf = zfun.filt_hanning(u, n=5)

df['uf'] = pd.Series(dict(zip(t,uf)))

# values in the array t are of type numpy.datetime64
# and you can convert to a regular datetime using
# pd.Timestamp(t[0]) but only for single values
#
# this gives and array of days since the first time
#a = ((t - t[0])/86400e9).astype(int)
#
# and this creates one of these objects
#b = np.datetime64(datetime(2012,1,1))
#
#tt = mdates.date2num(t)


#%% plotting

plt.close('all')

df.plot(subplots=True, grid=True)

df.ix[:,['u3d', 'uf']].plot(grid=True)


