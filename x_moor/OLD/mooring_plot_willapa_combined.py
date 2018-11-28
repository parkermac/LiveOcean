"""
Plots combined Willapa mooring records.
Focused just on surface water.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

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
indir = Ldir['LOo'] + 'extract/'
# choose the mooring extraction to plot

import willapa_fun as wfun
from importlib import reload
reload(wfun)

V_dict = dict()

for ww in ['0','1','2','3']:
    moor_file = ('moor_cas4_v1_lo6biom_willapa' + ww
        + '_hourly_2017.01.01_2017.12.31.nc')
    fn = indir + moor_file
    ds = nc.Dataset(fn)
    # create a dict into which to load everything
    V = dict()
    # gather some fields
    vn_list = ['salt', 'temp',
        'NO3', 'phytoplankton',
        'TIC', 'alkalinity',
        'PH', 'ARAG', 'z_rho']
    # gather data (including ocean_time and units
    for vn in vn_list:
        V[vn] = ds[vn][-1,:]
    V['svstr'] = ds['svstr'][:]
    V['ocean_time'] = ds['ocean_time'][:]
    V['lon'] = ds['lon_rho'][:]
    V['lat'] = ds['lat_rho'][:]
    ds.close()
    V_dict[ww] = V
    
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# load river files
riv_list = ['columbia', 'willapa', 'naselle']
riv_dict = dict()
for riv in riv_list:
    riv_dict[riv] = pd.read_pickle(Ldir['parent'] +
        'ptools_output/river/coastal_2017/' + riv + '.p')

#%% plotting
plt.close('all')
fs1=16 # fontsize for labels
fs2 = 12
lw=2

cdict = {'0':'b', '1':'g', '2':'orange', '3':'r'}

fig = plt.figure(figsize=(20,8))

V = V_dict['0']
mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
mdt = mdates.num2date(mdays) # list of datetimes of data
svstr = V['svstr'] # just use wund from the offshore mooring

nn_dict = {'salt': 1,
        'NO3': 2,
        'PH': 3,
        'temp': 5,
        'phytoplankton': 6,
        'ARAG': 7,
        'TIC': 4,
        'alkalinity': 8}

for vn in nn_dict:
    ax = fig.add_subplot(3,4,nn_dict[vn])
    
    for ww in ['0','1','2','3']:
        fld = V_dict[ww][vn]
        #ax.plot(mdt, wfun.fac_dict[vn] * zfun.filt_godin(fld), '-', color=cdict[ww])
        ax.plot(mdt, wfun.fac_dict[vn] * zfun.filt_hanning(fld, n=10*24),
            '-', color=cdict[ww], linewidth=lw)
    ax.grid(True)
    ax.set_xlim(mdt[0], mdt[-1])
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.text(.05, .85, wfun.tstr_dict[vn] + ' ' + wfun.units_dict[vn],
        horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.set_xticklabels([])
    ax.set_ylim(wfun.vlims_dict[vn])
    # # add depth labels
    # if vn=='salt':
    #     ax.text(.95, .05, ('Z = %d (m)' % (int(zbot))),
    #         horizontalalignment='right', color='b',
    #         transform=ax.transAxes, fontsize=fs2, fontweight='bold')
    #     ax.text(.95, .15, ('Z = %d (m)' % (int(zmid))),
    #         horizontalalignment='right', color='g',
    #         transform=ax.transAxes, fontsize=fs2, fontweight='bold')
    #     ax.text(.95, .25, ('Z = %d (m)' % (int(ztop))),
    #         horizontalalignment='right', color='r',
    #         transform=ax.transAxes, fontsize=fs2, fontweight='bold')

vn = 'svstr'
ax = fig.add_subplot(3,4,9)
fld = zfun.filt_AB8d(svstr)
ax.plot(mdt, fld, '-k', linewidth=lw)
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
ax.fill_between(mdt, fld, where=(fld<=0), color='slateblue', alpha=0.5)
ax.fill_between(mdt, fld, where=(fld>=0), color='lightsalmon', alpha=0.5)
ax.grid(True)
ax.set_xlim(mdt[0], mdt[-1])
ax.ticklabel_format(useOffset=False, axis='y')
ax.text(.05, .85, 'N-S Windstress, 8 day filter $(Pa)$',
    horizontalalignment='left', transform=ax.transAxes, fontsize=fs1)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_xticklabels([])
ax.set_ylim(-.1, .3)
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
ax.set_xlim(mdt[0], mdt[-1])
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.set_ylim(0, 25)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.set_xlabel('Date 2017')
ax.xaxis.set_tick_params(labelrotation=45)
ax.text(.05, .85, 'Columbia River Flow $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr, fontsize=fs1)

# rivers 2
ax = fig.add_subplot(3,4,11)
#
rr = riv_dict['willapa']
clr1 = 'teal'
rr = rr/1000
rr.plot(style='-', color=clr1, linewidth=lw)
ax.text(.05, .85, 'Willapa River Flow $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr1, fontsize=fs1)
#
rr = riv_dict['naselle']
clr2 = 'darkorange'
rr = rr/1000
rr.plot(style='-', color=clr2, linewidth=lw)
ax.text(.05, .7, 'Naselle River Flow $(1000\ m^{3}s^{-1})$',
    horizontalalignment='left', transform=ax.transAxes,
    color=clr2, fontsize=fs1)
#
ax.grid(True)
ax.set_xlim(mdt[0], mdt[-1])
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
for ww in ['0','1','2','3']:
    lon = V_dict[ww]['lon']
    lat = V_dict[ww]['lat']
    ax.plot(lon, lat, '*', color=cdict[ww], markersize=18)
ax.set_title('Stations')

fig.tight_layout()

if True:
    out_dir = Ldir['LOo']+'willapa_moorings/'
    Lfun.make_dir(out_dir)
    plt.savefig(out_dir + 'combined.png')
else:
    plt.show()

