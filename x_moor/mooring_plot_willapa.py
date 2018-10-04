"""
Plots mooring records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta

import pandas as pd

import os
import sys

pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun

pth = os.path.abspath('../plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'
# choose the mooring extraction to plot

import willapa_fun as wfun
from importlib import reload
reload(wfun)


nlay = 0 # which layer to choose

fs1=14 # fontsize for labels
fs2 = 12
lw=2

# load model mooring data
moor_file = 'cas4_v2_lo6biom_2017.01.01_2017.12.31/wbc_hourly.nc'
fn = indir + moor_file
ds = nc.Dataset(fn)
V = dict()
Vu = dict()
vn_list2 = ['svstr']
vn_list3 = ['salt', 'temp',
    'NO3', 'phytoplankton',
    'TIC', 'alkalinity',
    'PH', 'ARAG', 'z_rho']
vlist = vn_list2 + vn_list3
ocean_time = ds['ocean_time'][:]
for vn in vn_list2:
    V[vn] = ds[vn][:]
    try:
        Vu[vn] = ds[vn].units
    except AttributeError:
        Vu[vn] = ''
for vn in vn_list3:
    V[vn] = ds[vn][:, nlay].squeeze()
    try:
        Vu[vn] = ds[vn].units
    except AttributeError:
        Vu[vn] = ''
lon = ds['lon_rho'][:]
lat = ds['lat_rho'][:]
z = V['z_rho'].mean()
ds.close()
mdays = Lfun.modtime_to_mdate_vec(ocean_time)
mdt = mdates.num2date(mdays) # list of datetimes of data

# load river files
riv_list = ['columbia', 'willapa', 'naselle']
riv_dict = dict()
for riv in riv_list:
    riv_dict[riv] = pd.read_pickle(Ldir['parent'] +
        'ptools_output/river/coastal_2017/' + riv + '.p')
        
# observational data
dir0 = Ldir['parent'] + 'ptools_data/willapa/'
fn = dir0 + '2017Bay Center pCO2AverageBurke.xlsx'
df = pd.read_excel(fn)
# the original file timestamp does not have the year, and so
# pandas reads it as 1900. Here we try to make the correct time
dt = datetime(2017,1,1) - datetime(1900,1,1)
df['Date'] = df['Julian Day'] + dt# + timedelta(seconds=8*3600) # convert to UTC?
df = df.set_index('Date')
df.pop('Julian Day')
df = df.rename(columns={'Ω Arag':'Omega'})


#%% plotting
plt.close('all')
fig = plt.figure(figsize=(14,8))

nn_dict = {'salt': 1,
        'NO3': 2,
        'PH': 3,
        'temp': 5,
        'phytoplankton': 6,
        'ARAG': 7,
        'TIC': 4,
        'alkalinity': 8}
        
mdt0 = mdt[0]
mdt1 = mdt[-1]

for vn in nn_dict:
    ax = fig.add_subplot(3,4,nn_dict[vn])
    ax.plot(mdt, wfun.fac_dict[vn] * V[vn], '-b', linewidth=1)
    alph = .5
    if vn == 'ARAG':
        df.plot(y='Omega', ax=ax, legend=False, alpha=alph, color='r')
    elif vn == 'PH':
        df.plot(y='pH', ax=ax, legend=False, alpha=alph, color='r')
    elif vn == 'temp':
        df.plot(y='T', ax=ax, legend=False, alpha=alph, color='r')
    elif vn == 'salt':
        df.plot(y='S', ax=ax, legend=False, alpha=alph, color='r')
    ax.grid(True)
    ax.set_xlim(mdt0, mdt1)
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, wfun.tstr_dict[vn],
        horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
    ax.text(.05, .7, wfun.units_dict[vn],
        horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.set_xticklabels([])
    ax.set_ylim(wfun.vlims_dict[vn])
    # add depth labels
    if vn=='salt':
        ax.text(.95, .05, ('Model Z = %d (m)' % (int(z)) ),
            horizontalalignment='right', color='b',
            transform=ax.transAxes, fontsize=fs2, fontweight='bold')
        ax.text(.95, .15, 'Observed',
            horizontalalignment='right', color='r',
            transform=ax.transAxes, fontsize=fs2, fontweight='bold', alpha=alph)

vn = 'svstr'
ax = fig.add_subplot(3,4,9)
fld = zfun.filt_AB8d(V[vn][:])
ax.plot(mdt, fld, '-k')
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
ax.fill_between(mdt, fld, where=(fld<=0), color='slateblue', alpha=0.5)
ax.fill_between(mdt, fld, where=(fld>=0), color='lightsalmon', alpha=0.5)
ax.grid(True)
ax.set_xlim(mdt0, mdt1)
ax.ticklabel_format(useOffset=False, axis='y')
ax.text(.05, .85, 'N-S Windstress (8 day filter)',
    horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
ax.text(.05, .7, '$(Pa)$',
    horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_xticklabels([])
ax.set_ylim(-.05, .1)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.set_xlabel('Date 2017')
ax.xaxis.set_tick_params(labelrotation=45)

# rivers 1
ax = fig.add_subplot(3,4,10)
rr = riv_dict['columbia']
clr = 'darkorchid'
rr = rr/1000
rr.plot(style='-', color=clr, linewidth=lw)
ax.grid(True)
ax.set_xlim(mdt0, mdt1)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_ylim(0, 25)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.set_xlabel('Date 2017')
ax.xaxis.set_tick_params(labelrotation=45)
ax.text(.05, .85, 'Columbia River $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr, fontsize=fs1)

# rivers 2
ax = fig.add_subplot(3,4,11)
#
rr = riv_dict['willapa']
clr1 = 'teal'
rr = rr/1000
rr.plot(style='-', color=clr1, linewidth=lw)
ax.text(.05, .85, 'Willapa River $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr1, fontsize=fs1)
#
rr = riv_dict['naselle']
clr2 = 'darkorange'
rr = rr/1000
rr.plot(style='-', color=clr2, linewidth=lw)
ax.text(.05, .7, 'Naselle River $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr2, fontsize=fs1)
#
ax.grid(True)
ax.set_xlim(mdt0, mdt1)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_ylim(0, .5)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.set_xlabel('Date 2017')
ax.xaxis.set_tick_params(labelrotation=45)

# map
ax = fig.add_subplot(3,4,12)
pfun.add_coast(ax)
ax.axis([-124.6, -123.6, 46.3, 46.8])
pfun.dar(ax)
ax.plot(lon, lat, '*r')
sta = moor_file.split('/')[-1]
sta = sta.replace('.nc','')
ax.set_title(sta)

fig.tight_layout()
plt.show()


