"""
Plot bulk fluxes as a time series.  Similar to bulk plot, but cleaner,
designed for publication, and maybe focused only on one section.

We drop several panels that were in bulk_plot.py because we know now
what we need to about tidal energy flux, ssh, and net transport from
other plotting code.

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

import tef_fun
import flux_fun

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

testing = False

year = 2017
year_str = str(year)
indir0 = Ldir['LOo'] + 'tef/'
indir = indir0 + 'cas6_v3_lo8b_' + year_str + '.01.01_' + year_str + '.12.31/'

outdir = indir + 'bulk_plots_clean/'
Lfun.make_dir(outdir)

if testing:
    sect_list = ['jdf3']
else:
    sect_list = list(sect_df.index)

# PLOTTING
lw=3
fs=16
ms = 20
alpha = .5
qscl = 50
plt.rc('font', size=fs)
plt.close('all')

for sect_name in sect_list:
    
    fn = indir + 'bulk/' + sect_name + '.p'

    # some information about direction
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        if landward == 1:
            dir_str = 'Eastward'
        elif landward == -1:
            dir_str = 'Westward'
        a = [y0, y1]; a.sort()
        y0 = a[0]; y1 = a[1]
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        if landward == 1:
            dir_str = 'Northward'
        elif landward == -1:
            dir_str = 'Southward'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]

    # from warnings import filterwarnings
    # filterwarnings('ignore') # skip some warning messages

    bulk = pickle.load(open(fn, 'rb'))

    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']
    qnet = bulk['qnet_lp'] # net volume transport
    fnet = bulk['fnet_lp'] # net tidal energy flux
    ssh = bulk['ssh_lp'] # average SSH across the section, low-passed
    NT, NS = SS.shape

    # make vector and array times in days from start of the year
    dt = []
    for tt in ot:
        dt.append(Lfun.modtime_to_datetime(tt))
    td = []
    for tt in dt:
        ttt = tt - datetime(year,1,1)
        td.append(ttt.days + ttt.seconds/86400)
    td = np.array(td) # time in days from start of the year
    Time = td.reshape((NT,1)) * np.ones((1,NS)) # matrix version

    # separate out positive and negative transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan


    # ---------------------------------------------------------

    tef_df = flux_fun.get_fluxes(indir, sect_name, in_sign=1)
    tef_mean_df = tef_df.resample('1M').mean()
    # the above puts timestamps at the end of the month
    # so here we set it to the middle of each month because it is more
    # consistent with the averaging
    tef_mean_df.index -= timedelta(days=15)
    tef_mean_df.loc[:,'yd'] = tef_mean_df.index.dayofyear

    # ---------------------------------------------------------

    fig = plt.figure(figsize=(12,9))

    qlim_p = np.around(np.nanmax(1.3*QQp)/1000, 0)
    qlim_m = np.around(np.nanmax(-1.3*QQm)/1000, 0)
    qlim = np.max([qlim_p, qlim_m])

    # Salinity vs. Time (color by Transport)
    ax = fig.add_subplot(211)
    if False:
        ax.plot(Time[QQ>0], SS[QQ>0], 'or', alpha=alpha)
        ax.plot(Time[QQ<0], SS[QQ<0], 'ob', alpha=alpha)
    else:
        ax.scatter(Time, SS, s=np.abs(qscl*QQp/(1000*qlim)), c='r', alpha=alpha)
        ax.scatter(Time, SS, s=np.abs(qscl*QQm/(1000*qlim)), c='b', alpha=alpha)
    # add two-layer versions
    tef_mean_df.plot(x='yd', y = 'Sin', style='-*r', ax=ax, legend=False, mfc='w', ms=ms, lw=lw)
    tef_mean_df.plot(x='yd', y = 'Sout', style='-*b', ax=ax, legend=False, mfc='w', ms=ms, lw=lw)

    ax.set_title('Section = ' + sect_name + ': Positive is ' + dir_str)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Salinity')

    ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    # Tranport vs. Time
    ax = fig.add_subplot(212)
    if False:
        ax.plot(Time, QQp/1e3, 'or', alpha=alpha)
        ax.plot(Time, QQm/1e3, 'ob', alpha=alpha)
    else:
        ax.scatter(Time, QQp/1e3, s=np.abs(qscl*QQp/(1000*qlim)), c='r', alpha=alpha)
        ax.scatter(Time, QQm/1e3, s=np.abs(qscl*QQm/(1000*qlim)), c='b', alpha=alpha)
    # add two-layer versions
    this_yd = tef_mean_df.loc[:,'yd'].to_numpy()
    this_qin = tef_mean_df.loc[:,'Qin'].to_numpy()/1e3
    this_qout = tef_mean_df.loc[:,'Qout'].to_numpy()/1e3
    ax.plot(this_yd, this_qin, '-*r', mfc='w', ms=ms, lw=lw)
    ax.plot(this_yd, this_qout, '-*b', mfc='w', ms=ms, lw=lw)

    ax.plot([0, 366], [0,0], '-k', )

    ax.set_xlim(0,366)

    ax.set_ylim(-qlim, qlim)

    ax.grid(True)
    ax.set_xlabel('Yearday ' + year_str)
    ax.set_ylabel('Transport $[1000\ m^{3}s^{-1}]$')

    ax.text(.03, .95, '(b)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    if tef_mean_df['Sin'].mean() > tef_mean_df['Sout'].mean():
        pass
    else:
        print('Warning: sign logic breakdown! ' + sect_name)
    
    ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', color='r',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', color='b',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        
    fig.tight_layout()
    
    plt.savefig(outdir + sect_name + '.png')
    
    if testing:
        plt.show()
    else:
        plt.close()

plt.rcdefaults()
