#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:57:18 2018

@author: pm7

Code to plot a combination of tide records from ROMS
and PS_Tides.

"""

# setup
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from datetime import datetime, timedelta

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

def dn_to_dt(dn):
    # convert MATLAB datenum to datetime
    # there is a tiny roundoff error involved, amounts to 5e-6 sec
    dt = ( datetime.fromordinal(int(dn))
        + timedelta(days = np.mod(dn,1))
        - timedelta(days = 366) )
    return dt
    
def get_ps_tide(fn, offset1=0, offset2=0):
    # get data from a PS Tides mat file
    ps = scipy.io.loadmat(fn)
    t_pst = ps['Time_datenum_PST'].flatten()
    # convert to UTC
    t_utc = t_pst + 8/24
    # convert to array of python datetimes
    tt = []
    for t in t_utc:
        ttt = (dn_to_dt(t)
            + timedelta(days=offset1) + timedelta(days=offset2))
        ttt = ttt.replace(microsecond=0)
        tt.append(ttt)    
    T = np.array(tt)
    #
    h = ps['Height'].flatten()
    H = h - h.mean()
    #
    U = ps['Current'].flatten()
    
    LatLon = ps['LatLon'].flatten()
    # packed like:
    # [  47.80708, -122.73561,   47.77786, -122.69415]
    
    return T, H, U, LatLon
    
def get_mod_tide(fn, offset1=0, offset2=0):
    # get data from a ROMS mooring extraction
    ds = nc.Dataset(fn)
    #
    ot = ds['ocean_time'][:]
    tt = []
    for t in ot:
        ttt = (Lfun.modtime_to_datetime(t)
            + timedelta(days=offset1) + timedelta(days=offset2))
        tt.append(ttt)   
    T = np.array(tt)
    #
    zeta = ds['zeta'][:]
    H = zeta - zeta.mean()
    
    U = ds['ubar'][:]
    V = ds['vbar'][:]
    #
    ds.close()
    
    return T, H, U, V
    
# PS_tides input files
indir_ps = Ldir['parent'] + 'ptools_data/usrs/'
fn = indir_ps + 'SEG508_170401.mat'
T17, H17, U17, LatLon = get_ps_tide(fn, offset1=365, offset2=-10)
fn = indir_ps + 'SEG508_180401.mat'
T18, H18 , U18, LatLon= get_ps_tide(fn)

# model input file
indir_mod = Ldir['LOo'] + 'moor/'
fn = indir_mod + 'cas3_v1_lo6m_HCsill_backfill_2017.04.01_2017.06.01.nc'
TM, HM, UM, VM = get_mod_tide(fn, offset1=365, offset2=-10)

# convert all times from UTC to PDT
T17pdt = T17 - timedelta(seconds=7*3600)
T18pdt = T18 - timedelta(seconds=7*3600)
TMpdt = TM - timedelta(seconds=7*3600)

# plotting
plt.close('all')
fig, axes = plt.subplots(nrows=2, ncols=1,
    squeeze=False, sharex=True, figsize=(12,7))

ax = axes[0,0]
ax.plot(T17pdt, H17, '-c', label='PS Tides 17 (time adjusted)')
ax.plot(T18pdt, H18, '-b', label='PS Tides 18')
ax.plot(TMpdt, HM, '-r', label='Model (time adjusted)')
ax.set_xlim(T18[0],TM[-1])
ax.set_ylabel('SSH (m)')
ax.grid(True)
ax.legend()

ax = axes[1,0]
ax.plot(T17pdt, U17, '-c')
ax.plot(T18pdt, U18, '-b')
# define along-channel direction
th = 45 * np.pi/180
cc = np.cos(th)
ss = np.sin(th)
UU = cc*UM + ss*VM
VV = cc*VM - ss*UM
ax.plot(TMpdt, -UU, '-r')
ax.set_xlim(T18[0],TM[-1])
ax.set_xlabel('Date (PDT)')
ax.set_ylabel('U (m/s) Flood-Positive')
ax.grid(True)


plt.show()




