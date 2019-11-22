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
# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

year = 2017

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir0)
indir0 = indir0 + item + '/'
indir = indir0 + 'bulk/'
sect_list_raw = os.listdir(indir)
sect_list_raw.sort()
sect_list = [item for item in sect_list_raw if ('.p' in item)]
print(20*'=' + ' Processed Sections ' + 20*'=')
print(*sect_list, sep=", ")
print(61*'=')
# select which sections to process
my_choice = input('-- Input section to plot (e.g. sog5, or Return to plot all): ')
if len(my_choice)==0:
    # full list
    save_fig = True
else: # single item
    if (my_choice + '.p') in sect_list:
        sect_list = [my_choice + '.p']
        save_fig = False
    else:
        print('That section is not available')
        sys.exit()
outdir = indir0 + 'bulk_plots/'
Lfun.make_dir(outdir)
            
#plt.close('all')

sect_list = [item for item in sect_list if item.replace('.p','') in sect_df.index]

for snp in sect_list:
    
    sn = snp.replace('.p','')

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
        td_list.append(datetime(2017,1,1,0,0,0) + timedelta(days=t))
    tef_df = pd.DataFrame(index=td_list, columns=['Qin','Qout','Sin','Sout'])
    tef_df.loc[:,'Qin']=Qin
    tef_df.loc[:,'Qout']=Qout
    tef_df.loc[:,'Sin']=Sin
    tef_df.loc[:,'Sout']=Sout
    tef_mean_df = tef_df.resample('1M').mean()
    # the above puts timestamps at the end of the month
    # so here we set it to the middle of each month becasue it is more
    # consistent with the averaging
    tef_mean_df.index -= timedelta(days=15)
    tef_mean_df.loc[:,'yd'] = tef_mean_df.index.dayofyear
    

    # PLOTTING
    fig = plt.figure(figsize=(21,9))
    
    alpha = .5

    # Salinity vs. Time (size and color by Transport)
    ax = fig.add_subplot(2,3,1)
    Qscale = np.nanmean(np.abs(QQ))
    qf = 25
    ax.scatter(Time, SS, s=qf*np.abs(QQp/Qscale), c='r', alpha=alpha)
    ax.scatter(Time, SS, s=qf*np.abs(QQm/Qscale), c='b', alpha=alpha)
    # add two-layer versions
    if False:
        ax.plot(td, Sin, '-k', td, Sout, '--k')
    else:
        tef_mean_df.plot(x='yd', y = 'Sin', style='-ok', ax=ax, legend=False)
        tef_mean_df.plot(x='yd', y = 'Sout', style='--ok', ax=ax, legend=False)
    ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes, fontweight='bold')
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_ylabel('Salinity')
    # legend
    ax.scatter(.95, .2, s=qf, c='r', transform=ax.transAxes, alpha=alpha)
    ax.scatter(.95, .1, s=qf, c='b', transform=ax.transAxes, alpha=alpha)
    ax.text(.94, .2, 'Positive Q %d (m3/s)' % int(Qscale), color='r', fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.text(.94, .1, 'Negative Q %d (m3/s)' % int(Qscale), color='b', fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.set_title(indir0.split('/')[-2])
    
    # # Tidal energy flux vs. Time as second y-axis
    # ax = ax.twinx()
    # ax.plot(td, fnet/1e9, '-g', linewidth=2)
    # ax.set_ylabel('Energy Flux (GW)', color='g', alpha=alpha)
    # ax.set_ylim(bottom=0)
    # ax.set_xlim(0,366)
    
    # Tranport vs. Time
    ax = fig.add_subplot(2,3,4)
    ax.scatter(Time, QQp/1e3, s=qf*np.abs(QQp/Qscale), c='r', alpha=alpha)
    ax.scatter(Time, -QQm/1e3, s=qf*np.abs(QQm/Qscale), c='b', alpha=alpha)
    # add two-layer versions
    if False:
        ax.plot(td, Qin/1e3, '-k', td, -Qout/1e3, '--k')
    else:
        this_yd = tef_mean_df.loc[:,'yd'].values
        this_qin = tef_mean_df.loc[:,'Qin'].values/1e3
        this_qout = -tef_mean_df.loc[:,'Qout'].values/1e3
        # tef_mean_df.plot(x='yd', y = 'Qin', style='-ok', ax=ax, legend=False)
        # tef_mean_df.plot(x='yd', y = 'Qout', style='--ok', ax=ax, legend=False)
        ax.plot(this_yd, this_qin, '-ok')
        ax.plot(this_yd, this_qout, '--ok')
    ax.set_xlim(0,366)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(year))
    ax.set_ylabel('|Q| 1000 m3/s')
    
    # Tidal energy flux vs. Time as second y-axis
    ax = fig.add_subplot(3,3,2)
    ax.plot(td, fnet/1e9, '-g', linewidth=2)
    ax.set_ylabel('Energy Flux (GW)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(0,366)
    
    # Surface height
    ax = fig.add_subplot(3,3,5)
    ax.plot(td, ssh, '-b', linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_ylabel('SSH (m)')
    
    # Volume flux
    ax = fig.add_subplot(3,3,8)
    ax.plot(td, qnet/1e3, '-c', linewidth=2)
    ax.plot(td, Qnet/1e3, '--r', linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(year))
    ax.set_ylabel('Qnet 1000 m3/s')
    
    # Section location map
    ax = fig.add_subplot(1,3,3)
    ax.plot([x0, x1], [y0, y1], '-m', linewidth=3)
    ax.set_title(sn)
    pfun.add_coast(ax)
    pfun.dar(ax)
    aa = [x0-.7, x1+.7, y0-.5, y1+.5]
    ax.axis(aa)
    

    if save_fig:
        plt.savefig(outdir + sn + '.png')
        plt.close()
    else:
        plt.show()
