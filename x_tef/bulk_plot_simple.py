"""
Plot bulk fluxes as a time series.  This is like bulk_plot.py but
makes a simpler plot, designed for making points at Ocean Sciences 2020.
"""
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

import numpy as np
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta

import tef_fun
# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

year = 2017
year_str = str(year)
item = 'cas6_v3_lo8b_2017.01.01_2017.12.31'

# choose input and organize output
Ldir = Lfun.Lstart()
indir00 = Ldir['LOo'] + 'tef/'
indir0 = indir00 + item + '/'
indir = indir0 + 'bulk/'

outdir = indir00 + 'bulk_plots_simple/'
Lfun.make_dir(outdir)

plt.close('all')
for sn in ['hc1', 'ai1', 'mb3', 'jdf2', 'sji1', 'ss1', 'wb1']:
    qnet_sign = -1

    snp = sn + '.p'
    outname = outdir + sn + '_' + year_str + '.png'

    bulk = pickle.load(open(indir + snp, 'rb'))
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
        #ttt = tt- datetime(dt[0].year,1,1)
        ttt = tt - datetime(year,1,1) # hardwire for 2016.12.15 start
        td.append(ttt.days + ttt.seconds/86400)
    td = np.array(td) # time in days from start of the year
    Time = td.reshape((NT,1)) * np.ones((1,NS)) # matrix version

    # some information about direction
    x0, x1, y0, y1, landward = sect_df.loc[sn,:]    
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
    
    # separate out positive and negative transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan

    # form two-layer versions of Q and S
    Qin = np.nansum(QQp, axis=1)
    QSin = np.nansum(QQp*SS, axis=1)
    Sin = QSin/Qin
    Qout = np.nansum(QQm, axis=1)
    QSout = np.nansum(QQm*SS, axis=1)
    Sout = QSout/Qout
    # and find net transport to compare with qnet (should be identical)
    Qnet = np.nansum(QQ, axis=1)
    # RESULT: it is identical

    # 2019.11.14 make monthly averages
    import pandas as pd
    td_list = []
    for t in td:
        td_list.append(datetime(year,1,1,0,0,0) + timedelta(days=t))
    tef_df = pd.DataFrame(index=td_list, columns=['Qin','Qout','Sin','Sout','Qnet'])
    tef_df.loc[:,'Qin']=Qin
    tef_df.loc[:,'Qout']=Qout
    tef_df.loc[:,'Sin']=Sin
    tef_df.loc[:,'Sout']=Sout
    tef_df.loc[:,'Qnet'] = Qnet
    tef_mean_df = tef_df.resample('1M').mean()
    # the above puts timestamps at the end of the month
    # so here we set it to the middle of each month becasue it is more
    # consistent with the averaging
    tef_mean_df.index -= timedelta(days=15)
    tef_mean_df.loc[:,'yd'] = tef_mean_df.index.dayofyear

    # PLOTTING
    fig = plt.figure(figsize=(13,8))
    fs = 18

    # Salinity vs. Time (size and color by Transport)
    ax = fig.add_subplot(321)
    tef_mean_df.plot(x='yd', y = 'Sin', style='-ok', ax=ax, legend=False)
    tef_mean_df.plot(x='yd', y = 'Sout', style='--ok', ax=ax, legend=False)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.text(.05, .85, '(a) Sin and Sout', size=fs, transform=ax.transAxes)
    ax.tick_params(labelsize=fs-2)
    ax.set_xlabel('')

    # Tranport vs. Time
    ax = fig.add_subplot(323)
    this_yd = tef_mean_df.loc[:,'yd'].values
    this_qin = tef_mean_df.loc[:,'Qin'].values/1e3
    this_qout = -tef_mean_df.loc[:,'Qout'].values/1e3
    ax.plot(this_yd, this_qin, '-ok')
    ax.plot(this_yd, this_qout, '--ok')
    ax.set_xlim(0,366)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.text(.05, .05, '(b) Qin and Qout (1000 m3/s)', size=fs, transform=ax.transAxes)
    ax.tick_params(labelsize=fs-2)

    # Volume flux
    ax = fig.add_subplot(325)
    this_yd = tef_mean_df.loc[:,'yd'].values
    this_qnet = tef_mean_df.loc[:,'Qnet'].to_numpy()
    ax.plot(this_yd, qnet_sign * this_qnet/1000, '-ok', linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(year), size=fs)
    ax.text(.05, .05, '(c) Qriver (1000 m3/s)', size=fs, transform=ax.transAxes)
    ax.tick_params(labelsize=fs-2)

    # Section location map
    ax = fig.add_subplot(122)
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title('Section Location: ' + sn, size=fs)
    pfun.add_coast(ax)
    pfun.dar(ax)
    aa = [x0-.7, x1+.7, y0-.5, y1+.5]
    ax.axis(aa)
    ax.tick_params(labelsize=fs-4)

    plt.savefig(outname)
    
plt.show()
