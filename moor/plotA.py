#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:13:28 2016

@author: PM5

Plots mooring records for analytical runs like aestus1.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
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

#%% load and organize data

v2_list = []
v3_list_rho = []
v3_list_w = []
ds = nc.Dataset(fn)
for vv in ds.variables:
    vdim = ds.variables[vv].dimensions
    if ( ('ocean_time' in vdim)
        and ('s_rho' not in vdim)
        and ('s_w' not in vdim)
        and (vv != 'ocean_time') ):
        v2_list.append(vv)
    elif ( ('ocean_time' in vdim) and ('s_rho' in vdim) ):
        v3_list_rho.append(vv)
    elif ( ('ocean_time' in vdim) and ('s_w' in vdim) ):
        v3_list_w.append(vv)

# load everything into a dict
V = dict()

list_to_plot = v3_list_rho + v3_list_w + v2_list
list_to_plot.remove('temp')
list_to_plot.remove('shflux')
list_to_plot.remove('sustr')
list_to_plot.remove('svstr')
list_to_plot.remove('hh')

for vv in list_to_plot:
    V[vv] = ds[vv][:]

V['ocean_time'] = ds['ocean_time'][:]

ds.close()

#%% plotting

NP = len(list_to_plot)

NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
NC = np.ceil(np.sqrt(NP)).astype(int)

if NR*NC - NP >= NR:
    NC = NC-1

fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
mdt = mdates.num2date(mdays) # list of datetimes of data

cc = 0
nmid = 20
for vn in list_to_plot:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    if V[vn].ndim == 2:
        ax.plot(days, V[vn][-1,:], '-r')
        ax.plot(days, V[vn][nmid,:],'-g')
        ax.plot(days, V[vn][0,:], '-b')
    elif V[vn].ndim == 1:
        ax.plot(days, V[vn])

    # general case
    ax.set_xlim(days[0], days[-1])
    if ir == NR-1:
        ax.set_xlabel('Days')
    aa = ax.get_ylim()
    ax.set_ylim(aa)
        
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, vn,
            horizontalalignment='left',
            transform=ax.transAxes, fontweight='bold', fontsize=24)
    cc += 1

plt.show()
