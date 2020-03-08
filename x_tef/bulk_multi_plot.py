"""
Plot bulk fluxes as a time series.

Plots two or more sections on the same axes.
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
item = 'cas6_v3_lo8b_'+year_str+'.01.01_'+year_str+'.12.31'

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tef/' + item + '/bulk/'

sect_list = ['ai3', 'wb1']
color_list = [('orangered','darkorange'),('mediumblue','darkturquoise')]

plt.close('all')
# PLOTTING
fig = plt.figure(figsize=(21,9))

ax1 = fig.add_subplot(2,3,1)
ax2 = fig.add_subplot(2,3,2)
ax3 = fig.add_subplot(2,3,4)
ax4 = fig.add_subplot(2,3,5)
ax4a = ax4.twinx()
axm = fig.add_subplot(1,3,3)

ii = 0
for sn in sect_list:
    
    incolor = color_list[ii][0]
    outcolor = color_list[ii][1]
    
    bulk = pickle.load(open(indir + sn + '.p', 'rb'))
    QQ = bulk['QQ']
    SS = bulk['SS']
    ot = bulk['ot']
    qnet = bulk['qnet_lp'] # net volume transport
    fnet = bulk['fnet_lp'] # net tidal energy flux
    ssh = bulk['ssh_lp'] # average SSH across the section, low-passed
    NT, NS = SS.shape
    
    # separate out positive and negative transports
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan

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
    
    alpha = .5

    # Salinity vs. Time (size and color by Transport)
    ax = ax1
    Qscale = 1e4
    qf = 25
    ax.scatter(Time, SS, s=qf*np.abs(QQp/Qscale), c=incolor, alpha=alpha)
    ax.scatter(Time, SS, s=qf*np.abs(QQm/Qscale), c=outcolor, alpha=alpha)
    ax.text(0.05, 0.1, 'Positive is ' + dir_str, transform=ax.transAxes, fontweight='bold')
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_ylabel('Salinity')
    # legend
    ax.scatter(.95, .2, s=qf, c=incolor, transform=ax.transAxes, alpha=alpha)
    ax.scatter(.95, .1, s=qf, c=outcolor, transform=ax.transAxes, alpha=alpha)
    ax.text(.94, .2, 'Positive Q %d (m3/s)' % int(Qscale), color=incolor, fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.text(.94, .1, 'Negative Q %d (m3/s)' % int(Qscale), color=outcolor, fontweight='bold',
        horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.set_title(item)
    
    # Surface height
    ax = ax2
    ax.plot(td, ssh, '-', color=incolor, linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_ylabel('SSH (m)')
    
    # Tranport vs. Time
    ax = ax3
    ax.scatter(Time, QQp/1e3, s=qf*np.abs(QQp/Qscale), c=incolor, alpha=alpha)
    ax.scatter(Time, -QQm/1e3, s=qf*np.abs(QQm/Qscale), c=outcolor, alpha=alpha)
    ax.set_xlim(0,366)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(year))
    ax.set_ylabel('|Q| 1000 m3/s')
    
    # Volume flux
    ax = ax4
    ax.plot(td, qnet/1e3, '-', color=incolor, linewidth=2)
    ax.set_xlim(0,366)
    ax.grid(True)
    ax.set_xlabel('Days from 1/1/' + str(year))
    ax.set_ylabel('Qnet 1000 m3/s')
    
    # Tidal energy flux vs. Time as second y-axis
    ax = ax4a
    ax.plot(td, fnet/1e9, '-', color=outcolor, linewidth=2)
    ax.set_ylabel('Energy Flux (GW)')
    #ax.set_ylim(bottom=0)
    ax.set_xlim(0,366)
    
    
    # Section location map
    ax = axm
    ax.plot([x0, x1], [y0, y1], '-', color=incolor, linewidth=3)
    ax.set_title(sn)
    pfun.add_coast(ax)
    pfun.dar(ax)
    aa = [x0-.7, x1+.7, y0-.5, y1+.5]
    ax.axis(aa)
    
    ii += 1
    

plt.show()
