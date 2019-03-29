"""
Plot bulk fluxes as a time series.
"""
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import numpy as np
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta

import tef_fun

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

indir0 = ('/Users/pm7/Documents/LiveOcean_output/tef/' +
            'cas5_v3_lo8_2017.01.01_2017.06.20/')
            
indir = indir0 + 'bulk/'
            
outdir = indir0 + 'bulk_plots/'

if True: # plot all bulk files
    Lfun.make_dir(outdir, clean=True)
    LList = [item for item in os.listdir(indir) if ('.p' in item)]
    save_fig = True
else: # override
    snp = Lfun.choose_item(indir, tag='.p')
    Lfun.make_dir(outdir)
    LList = [snp]
    save_fig = False
            
for snp in LList:
    
    sn = snp.replace('.p','')

    bulk = pickle.load(open(indir + snp, 'rb'))
    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']
    qnet = bulk['qnet_lp'] # net volume transport
    fnet = bulk['fnet_lp'] # net tidal energy flux
    NT, NS = SS.shape

    # make vector and array times in days from start of the year
    dt = []
    for tt in ot:
        dt.append(Lfun.modtime_to_datetime(tt))
    td = []
    for tt in dt:
        ttt = tt- datetime(dt[0].year,1,1)
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

    # PLOTTING
    plt.close('all')
    fig = plt.figure(figsize=(14,7))
    
    alpha = .5

    # Salinity vs. Time (size and color by Transport)
    ax = plt.subplot2grid((2,3), (0,0), colspan=2)
    Qscale = np.nanmean(np.abs(QQ))
    qf = 25
    ax.scatter(Time, SS, s=qf*np.abs(QQp/Qscale), c='r', alpha=alpha)
    ax.scatter(Time, SS, s=qf*np.abs(QQm/Qscale), c='b', alpha=alpha)
    # add two-layer versions
    ax.plot(td, Sin, '-k', td, Sout, '--k')
    ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes, fontweight='bold')
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_ylabel('Salinity')
    # legend
    ax.scatter(.95, .2, s=qf, c='r', transform=ax.transAxes, alpha=alpha)
    ax.scatter(.95, .1, s=qf, c='b', transform=ax.transAxes, alpha=alpha)
    ax.text(.94, .2, 'Positive Q %d $(m^{3}s^{-1})$' % int(Qscale), color='r', fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.text(.94, .1, 'Negative Q %d $(m^{3}s^{-1})$' % int(Qscale), color='b', fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    # Tidal energy flux vs. Time as second y-axis
    ax = ax.twinx()
    ax.plot(td, fnet/1e9, '-g', linewidth=2)
    ax.set_ylabel('Energy Flux (GW)', color='g', alpha=alpha)
    ax.set_ylim(bottom=0)
    ax.set_xlim(0,366)
    
    # Tranport vs. Time
    ax = plt.subplot2grid((2,3), (1,0), colspan=2)
    ax.scatter(Time, QQp/1e3, s=qf*np.abs(QQp/Qscale), c='r', alpha=alpha)
    ax.scatter(Time, -QQm/1e3, s=qf*np.abs(QQm/Qscale), c='b', alpha=alpha)
    # add two-layer versions
    ax.plot(td, Qin/1e3, '-k', td, -Qout/1e3, '--k')
    ax.set_xlim(0,366)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(dt[0].year))
    ax.set_ylabel('|Q| $(10^{3} m^{3}s^{-1})$')
    # Tidal energy flux vs. Time as second y-axis
    ax = ax.twinx()
    ax.plot(td, fnet/1e9, '-g', linewidth=2)
    ax.set_ylabel('Energy Flux (GW)', color='g', alpha=alpha)
    ax.set_ylim(bottom=0)
    ax.set_xlim(0,366)
    
    # Section location map
    ax = fig.add_subplot(244)
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sn)
    pfun.add_coast(ax)
    pfun.dar(ax)
    aa = [x0-.7, x1+.7, y0-.5, y1+.5]
    ax.axis(aa)
    
    # Test of volume flux consistency
    ax = fig.add_subplot(248)
    ax.plot(td, Qnet/1e3, '-c', td, qnet/1e3, '--m', linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(dt[0].year))
    ax.set_ylabel('Qnet $(10^{3} m^{3}s^{-1})$')


    if save_fig:
        plt.savefig(outdir + sn + '.png')
        plt.close()
    else:
        plt.show()
