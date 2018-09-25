"""
Plot the mean of many TEF extractions on a thalweg section.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir = Ldir['LOo'] + 'tef/'
if False:
    print('\nSelect an Extraction to plot:\n')
    List = os.listdir(indir)
    List.sort()
    NL = len(List)
    Ldict = dict(zip(range(NL), List))
    for ii in range(NL):
        print(str(ii) + ': ' + List[ii])
    my_ii = int(input('-- Input number: '))
    Litem = Ldict[my_ii]
else:
    Litem = 'cas4_v1_lo6biom_2017.01.01_2017.12.31'
print('\nProcessing ' + Litem + '\n')
Indir = indir + Litem + '/'

plt.close('all')
# plotting
lw=2
fig = plt.figure(figsize=(13,8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

cc = 0
for sect_list in [['jdf1','jdf2','jdf3','jdf4',
    'ai1', 'ai2', 'ai3','ai4',
    'mb1','mb2','mb3','mb4','mb5',
    'tn1','tn2','tn3',
    'ss1','ss2','ss3'],
    ['jdf1','jdf2','jdf3','jdf4',
    'sji1', 'sji2', 'sog1','sog2',
    'sog3','sog4']]:

    NS = len(sect_list)
    dd = np.nan * np.ones(NS)
    qin = np.nan * np.ones(NS)
    qout = np.nan * np.ones(NS)
    qsin = np.nan * np.ones(NS)
    qsout = np.nan * np.ones(NS)
    qin_abs = np.nan * np.ones(NS)
    qout_abs = np.nan * np.ones(NS)
    qsin_abs = np.nan * np.ones(NS)
    qsout_abs = np.nan * np.ones(NS)
    sin = np.nan * np.ones(NS)
    sout = np.nan * np.ones(NS)
    xs = np.nan * np.ones(NS)
    ys = np.nan * np.ones(NS)
    counter = 0
    for sect_name in sect_list:    
        print('** ' + sect_name + ' **')    
        x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
        lon = (x0+x1)/2
        lat = (y0+y1)/2    
        if counter == 0:
            lon0 = lon
            lat0 = lat        
        xs[counter], ys[counter] = zfun.ll2xy(lon, lat, lon0, lat0)
        fn = Indir + sect_name + '.p'
        Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(fn)        
        qin[counter] = np.nanmean(Qi[:,0]/1e3)
        qout[counter] = np.nanmean(Qi[:,1]/1e3)        
        qsin[counter] = np.nanmean(Fi[:,0]/1e3)
        qsout[counter] = np.nanmean(Fi[:,1]/1e3)    
        qin_abs[counter] = np.nanmean(np.abs(Qi[:,0]))
        qout_abs[counter] = np.nanmean(np.abs(Qi[:,1]))        
        qsin_abs[counter] = np.nanmean(np.abs(Fi[:,0]))
        qsout_abs[counter] = np.nanmean(np.abs(Fi[:,1]))    
        sin[counter] = qsin_abs[counter]/qin_abs[counter]
        sout[counter] = qsout_abs[counter]/qout_abs[counter]        
        counter += 1
    # create a distance vector
    dx = np.diff(xs)
    dy = np.diff(ys)
    dd = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(NS)
    dist[1:] = np.cumsum(dd/1000)
    
    if cc == 0:
        lstr = 'JdF to South Sound'
        lcol = 'olive'
    elif cc == 1:
        lstr = 'JdF to Strait of Georgia'
        lcol = 'darkorange'
        
    ax1.plot(dist,qin,'-o', color=lcol,linewidth=lw, label=lstr)
    ax1.set_xlim(0,400)
    ax1.grid(True)
    ax1.set_ylabel('Qin [1000 m3/s]')
    ax1.legend()
    counter = 0
    for sn in sect_list:
        sn = sn.upper()
        ax1.text(dist[counter], qin[counter]+20, sn, rotation=45)
        counter += 1
        
    ax2.plot(dist,sin,'-+', color=lcol ,linewidth=lw, label='Sin ' + lstr)
    ax2.plot(dist,sout,'-+', color=lcol,linewidth=lw, label='Sout ' + lstr, alpha=.5)
    ax2.set_xlim(0,400)
    ax2.grid(True)
    ax2.set_xlabel('Distance from Mouth (km)')
    ax2.set_ylabel('Sin and Sout [psu]')
    ax2.legend()
    
    cc += 1

plt.show()
    