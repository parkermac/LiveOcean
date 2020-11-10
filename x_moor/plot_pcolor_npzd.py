"""
Plots mooring records as pcolormesh fields C(t,z).

Currently oriented around NPZD variables.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from datetime import datetime

import os; import sys
sys.path.append(os.path.abspath('../alpha'))

import Lfun
import zfun
low_pass = True

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

# choose the mooring extraction to plot
if False:
    item = Lfun.choose_item(indir0)
    indir = indir0 + item + '/'
    infile = Lfun.choose_item(indir, tag='.nc')
else:
    indir = (Ldir['LOo'] + 'moor/'
        +'cas6_v3_lo8b_2017.01.01_2019.12.31/')
    infile = 'HCB003_hourly.nc'
fn = indir + infile
# NOTE: these are packed (t,z)

ds = nc.Dataset(fn)
# time
ot_vec = ds['ocean_time'][:].data
dt_list = []
for ot in ot_vec:
    dt_list.append(Lfun.modtime_to_datetime(ot))
mdays = mdates.date2num(dt_list)
days = mdates.num2date(mdays)
Days = days[12::24]
# space
zr =  ds['z_rho'][:].data
zw = ds['z_w'][:].data
Zr = zr.mean(axis=0)
Zw = zw.mean(axis=0)

# colormaps
cmap_dict = {'temp': 'RdYlBu_r',
        'salt': 'Spectral_r',
        'phytoplankton': 'ocean_r',
        'zooplankton': 'cubehelix_r',
        'detritus': 'cividis_r',
        'Ldetritus': 'copper_r',
        'oxygen': 'rainbow_r',
        'TIC': 'cool',
        'alkalinity': 'cool',
        'NO3': 'cool',
        'Ntotal': 'cool'}

if True:
    vn_list = ['phytoplankton', 'zooplankton',
        'detritus', 'Ldetritus', 'NO3', 'Ntotal']
else:
    vn_list = ['salt', 'temp']
NR = len(vn_list)

plt.close('all')
fig = plt.figure(figsize=(13,8))
fs = 14

rr = 1
for vn in vn_list:
    if vn == 'Ntotal':
        Nvn_list = ['phytoplankton', 'zooplankton', 'NO3',
            'detritus', 'Ldetritus']
        ii = 0
        for Nvn in Nvn_list:
            if ii == 0:
                v = ds[Nvn][:].data
            else:
                v = v + ds[Nvn][:].data
            ii += 1
        Ntotal = v.copy()
        units = 'millimole_nitrogen meter-3'
    else:
        v = ds[vn][:].data
        if vn == 'salt':
            units = 'g/kg'
        else:
            units = ds[vn].units
    V = zfun.filt_godin_mat(v)[12::24,:].T
    cmap = cmap_dict[vn]
    ax = fig.add_subplot(NR, 1, rr)
    cs = ax.pcolormesh(Days, Zr, V, cmap=cmap)
    if vn == 'NO3':
        Nlow, Nhi = cs.get_clim()
    if vn == 'Ntotal':
        cs.set_clim(vmin=Nlow, vmax=Nhi)
    plt.vlines([Days[365],Days[2*365]],Zr[0],Zr[-1])
    if True:
        cb = fig.colorbar(cs, aspect=20/NR)
    else:
        # Inset colorbar
        cbaxes = inset_axes(ax, width="25%", height="8%", loc='lower right', borderpad=3)
        cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
        cb.ax.tick_params(labelsize=fs-2)
    ax.text(.05, .1, vn + ': (' + units.replace('_',' ') + ')',
        transform=ax.transAxes,
        style='italic', size=fs+1)
    if rr == 1:
        ax.set_title(infile, size=fs+2)
    if rr < NR:
        ax.set_xticklabels([])
    if rr == NR:
        ax.set_xlabel('Date', size=fs)
    ax.set_ylabel('Z (m)', size=fs)
    ax.tick_params(labelsize=fs-2) # tick labels
    rr += 1

if ('Ntotal' in vn_list) and False:
    # plot total watercolumn N
    dz = np.diff(zw, axis=1)
    Ntotal_A = np.sum(Ntotal*dz, axis=1)
    Ntotal_A_lp = zfun.filt_godin(Ntotal_A)
    Ntotal_A_lp = Ntotal_A_lp[12::24]
    fig2 = plt.figure(figsize=(13,8))
    ax = fig2.add_subplot(111)
    ax.plot(Days, Ntotal_A_lp, '-k')
    ax.set_xlim((Days[0], Days[-1]))
    ax.set_xlabel('Date', size=fs)
    ax.text(.05, .05, 'Integrated Total N: millimole nitrogen meter-3',
        transform=ax.transAxes,
        style='italic', size=fs)
    ax.set_title(infile, size=fs+2)
    ax.tick_params(labelsize=fs-2) # tick labels

plt.show()

