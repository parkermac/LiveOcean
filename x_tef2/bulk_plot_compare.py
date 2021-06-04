"""
Plot bulk fluxes as a time series.  Meant to focus on side-by-side
comparisons of two variations (e.g. one with smaller tides)

"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart()

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from datetime import datetime, timedelta
import netCDF4 as nc

import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

sect_name = 'ai1'
testing = True

indir0 = Ldir['LOo'] + 'tef2/'
indir_a = indir0 + 'cas6_v3_lo8b_2018.01.01_2018.12.31/'
indir_b = indir0 + 'cas6_v3t075_lo8_2018.07.01_2018.08.31/'

outdir = indir0 + 'bulk_plots_compare/'
Lfun.make_dir(outdir)

# PLOTTING
lw=2
fs=16
ms = 20
alpha = .2
qscl = 50
plt.rc('font', size=fs)
plt.close('all')

fig = plt.figure(figsize=(18,11))

ii = 1
for indir in [indir_a, indir_b]:
        

    # ---------------------------------------------------------

    tef_df, in_sign = flux_fun.get_fluxes(indir, sect_name)
    tef_df['Qout'] = -tef_df['Qout']/1000
    tef_df['Qin'] = tef_df['Qin']/1000
    tef_df['DS'] = tef_df['Sin'] - tef_df['Sout']
    tef_df['Qtide'] = tef_df['Qtide']/1000
    
    # some information about direction
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        if in_sign == 1:
            dir_str = 'Eastward'
        elif in_sign == -1:
            dir_str = 'Westward'
        a = [y0, y1]; a.sort()
        y0 = a[0]; y1 = a[1]
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        if in_sign == 1:
            dir_str = 'Northward'
        elif in_sign == -1:
            dir_str = 'Southward'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]
    
    # ---------------------------------------------------------
    
    dt0 = datetime(2018,7,1)
    dt1 = datetime(2018,8,31)
        
    # Salinity vs. Time
    ax = fig.add_subplot(3,2,ii)
    tef_df['Sin'].plot(c='r', lw=lw, ax=ax)
    tef_df['Sout'].plot(c='b', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(30,33)
    ax.set_title(indir.split('/')[-2])
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Salinity')
    ax.text(.03, .95, '(a) Section = ' + sect_name, va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .95, '$S_{in}$', ha='right', va='top', weight='bold', color='r',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .05, '$S_{out}$', ha='right', va='bottom', weight='bold', color='b',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    # Tidal Transport vs. Time
    ax = fig.add_subplot(3,2,ii+2)
    tef_df['Qtide'].plot(c='k', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(0,350)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.text(.03, .95, '(b) $Q_{tide}\ [10^{3}m^{3}s^{-1}]$',
        va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    # Tranport vs. Time
    ax = fig.add_subplot(3,2,ii+4)
    tef_df['Qin'].plot(c='r', lw=lw, ax=ax)
    tef_df['Qout'].plot(c='b', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(0,70)
    ax.grid(True)
    ax.set_xlabel('Date ' + str(dt0.year))
    ax.set_ylabel('')
    ax.text(.03, .95, '(c) $Q_{in}\ Q_{out}\ [10^{3}m^{3}s^{-1}]$',
        va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
    ii += 1
    
fig.tight_layout()
plt.savefig(outdir + sect_name + '.png')
plt.show()

plt.rcdefaults()
