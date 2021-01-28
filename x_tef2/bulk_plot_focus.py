"""
Plot bulk fluxes as a time series.  Meant to focus on making
plots for my V-PECS talk and an NSF proposal.

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

testing = True
if testing:
    sect_list = ['ai1', 'ai4']
else:
    sect_list = list(sect_df.index)

year = 2018
year_str = str(year)
indir0 = Ldir['LOo'] + 'tef2/'
indir = indir0 + 'cas6_v3_lo8b_' + year_str + '.01.01_' + year_str + '.12.31/'

outdir = indir + 'bulk_plots_focus/'
Lfun.make_dir(outdir)

# PLOTTING
lw=2
fs=16
ms = 20
alpha = .2
qscl = 50
plt.rc('font', size=fs)
plt.close('all')

def spring_neap(ax, txt=False):
    ax.axvspan(dt0_spring, dt1_spring, alpha=alpha, color='orange')
    ax.axvspan(dt0_neap, dt1_neap, alpha=alpha, color='green')
    if txt:
        ax.text(dt0_spring+timedelta(days=3.5), ax.get_ylim()[0]+.2,'SPRING',
            weight='bold', c='orange', ha='center', size=1.2*fs)
        ax.text(dt0_neap+timedelta(days=3.5), ax.get_ylim()[0]+.2,'NEAP',
            weight='bold', c='g', ha='center', size=1.2*fs)
        
    

for sect_name in sect_list:
        

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

    fig = plt.figure(figsize=(12,11))
    
    dt0 = datetime(2018,4,1)
    dt1 = datetime(2018,6,30)
    
    dt0_spring = datetime(2018,5,15)
    dt1_spring = dt0_spring + timedelta(days=7)
    dt0_neap = datetime(2018,6,5)
    dt1_neap = dt0_neap + timedelta(days=7)
    
    # Salinity vs. Time
    ax = fig.add_subplot(311)
    tef_df['Sin'].plot(c='r', lw=lw, ax=ax)
    tef_df['Sout'].plot(c='b', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_title('Section = ' + sect_name)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Salinity')
    ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    spring_neap(ax, txt=True)
    ax.text(.97, .95, '$S_{in}$', ha='right', va='top', weight='bold', color='r',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .05, '$S_{out}$', ha='right', va='bottom', weight='bold', color='b',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    # Tidal Transport vs. Time
    ax = fig.add_subplot(312)
    tef_df['Qtide'].plot(c='k', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.text(.03, .95, '(b) $Q_{tide}\ [10^{3}m^{3}s^{-1}]$',
        va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    spring_neap(ax)

    # Tranport vs. Time
    ax = fig.add_subplot(313)
    tef_df['Qin'].plot(c='r', lw=lw, ax=ax)
    tef_df['Qout'].plot(c='b', lw=lw, ax=ax)
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xlabel('Date ' + year_str)
    ax.set_ylabel('')
    ax.text(.03, .95, '(c) $Q_{in}\ Q_{out}\ [10^{3}m^{3}s^{-1}]$',
        va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    spring_neap(ax)
        
    fig.tight_layout()
    plt.savefig(outdir + sect_name + '.png')
    if testing:
        plt.show()
    else:
        plt.close()

plt.rcdefaults()
