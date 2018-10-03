"""
Plots mooring records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
from importlib import reload
reload(Lfun)

auto_lims = True

# set limits
lim_dict = {'temp': (0, 20),
        'salt': (25, 35),
        'NO3': (-1, 48),
        'phytoplankton': (-1, 14),
        'zooplankton': (-.1, 1),
        'detritus': (-.1, 2.5),
        'Ldetritus': (-.01, .15),
        'oxygen': (-10, 350),
        'TIC': (1800, 2600),
        'alkalinity': (1800, 2600)}

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

# choose the mooring extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
infile = Lfun.choose_item(indir, tag='.nc')
fn = indir + infile

#%% load and organize data

v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']
v1_list = ['ocean_time']
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

# create a dict into which to load everything
V = dict()
# and a dict of units
Vu = dict()

#choose what to plot
list_to_plot = v3_list_rho + v3_list_w + v2_list
#list_to_plot = v3_list_rho

# hand edit variables not to look at
for v in ['CaCO3']:#, 'PH', 'ARAG']:
    try:
        list_to_plot.remove(v)
    except ValueError:
        pass

ltp = list_to_plot.copy()
ltp.append('ocean_time')
# gather data (including ocean_time and units
for vv in ltp:
    V[vv] = ds[vv][:]
    try:
        Vu[vv] = ds[vv].units
    except AttributeError:
        Vu[vv] = ''

ds.close()

#%% plotting

plt.close('all')

NP = len(list_to_plot)
NR, NC = zfun.get_rc(NP)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(13,8),
                         squeeze=False, sharex=True)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

cc = 0
nmid = round(V['z_rho'].shape[1]/2)
N = V['z_rho'].shape[1]
nbot = 0
nmid = nmid
ntop = N-1

for vn in list_to_plot:
    
    ir, ic = zfun.get_irc(cc, NC)
    ax = axes[ir, ic]
    if True: # raw
        if V[vn].ndim == 2:
            ax.plot(days, V[vn][:, ntop], '-r')
            ax.plot(days, V[vn][:, nmid],'-g')
            ax.plot(days, V[vn][:, nbot], '-b')
        elif V[vn].ndim == 1:
            ax.plot(days, V[vn])
    else: # filtered (e.g. tidally_averaged)
        if V[vn].ndim == 2:
            ax.plot(days, zfun.filt_godin(V[vn][:, ntop]), '-r')
            ax.plot(days, zfun.filt_godin(V[vn][:, nmid]),'-g')
            ax.plot(days, zfun.filt_godin(V[vn][:, nbot]), '-b')
        elif V[vn].ndim == 1:
            ax.plot(days, zfun.filt_godin(V[vn]))
    
    try:
        if not auto_lims:
            ax.set_ylim(lim_dict[vn][0], lim_dict[vn][1])
    except KeyError:
        pass
        
    ax.grid(True)
    ax.set_xlim(days[0], days[-1])

    if ir == NR-1:
        ax.set_xlabel('Days')

    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, vn,
            horizontalalignment='left',
            transform=ax.transAxes)
    ax.text(.05, .75, Vu[vn],
            horizontalalignment='left',
            transform=ax.transAxes)
    cc += 1

plt.suptitle(fn)

plt.show()
