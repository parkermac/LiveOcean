# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:24:16 2016

@author: PM5
"""

"""
Plots mooring records.  Created for the Pacific Anomalies Workshop 2016_01
"""

# setup
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import numpy as np

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)

from importlib import reload
import Lfun
reload(Lfun)

import zfun  # plotting functions
import matfun  # functions for working with mat files

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the data to plot
if False:
    print('\n%s\n' % '** Choose mooring file to plot **')
    m_list_raw = os.listdir(indir)
    m_list = []
    for m in m_list_raw:
        if m[-2:] == '.p':
            m_list.append(m)
    Npt = len(m_list)
    m_dict = dict(zip(range(Npt), m_list))
    for npt in range(Npt):
        print(str(npt) + ': ' + m_list[npt])
    my_npt = int(input('-- Input number -- '))
    inname = m_dict[my_npt]
else:
    inname = 'cascadia1_base_lo1_47N124.5W_low_pass_2013.01.02_2015.12.31.p'
  
import pickle
V, v1_list, v2_list, v3_list, G, S = pickle.load( open( indir + inname, 'rb' ) )

import matplotlib.pyplot as plt
plt.close()
NR = 3; NC = 4
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])

year = 2013 + (mdays - mdates.date2num(datetime(2013,1,1)))/356

cc = 0
nmid = round(V['salt'].shape[0]/2)
nfilt = 20
for vn in v3_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(year, zfun.filt_hanning(V[vn][-1,:],n=nfilt), '-r')
    ax.plot(year, zfun.filt_hanning(V[vn][nmid,:],n=nfilt),'-g')
    ax.plot(year, zfun.filt_hanning(V[vn][0,:],n=nfilt), '-b')
    ax.set_xlim(2013, 2016)
    ax.ticklabel_format(useOffset=False, axis='x')
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_xticks([2013.5, 2014.5, 2015.5])
    ax.set_xticklabels([2013, 2014, 2015])
    ax.set_title(vn)
    yl = ax.get_ylim()
    ax.plot([2014, 2014],yl,'-k')
    ax.plot([2015, 2015],yl,'-k')
    ax.set_ylim(yl)
    cc += 1
for vn in v2_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(year, zfun.filt_hanning(V[vn],n=nfilt))
    ax.set_xlim(2013, 2016)
    ax.ticklabel_format(useOffset=False, axis='x')
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_xticks([2013.5, 2014.5, 2015.5])
    ax.set_xticklabels([2013, 2014, 2015])
    ax.set_title(vn)
    yl = ax.get_ylim()
    ax.plot([2014, 2014],yl,'-k')
    ax.plot([2015, 2015],yl,'-k')
    ax.set_ylim(yl)
    cc += 1
plt.show()
