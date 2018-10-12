"""
Plot the mean of many TEF extractions on a thalweg section.
"""

# setup
import numpy as np
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


indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'

# a structure to hold results for future use
ThalMean = dict()

cc = 0
for sect_list in [['jdf1','jdf2','jdf3','jdf4', # JdF to South Sound
            'ai1', 'ai2', 'ai3','ai4',
            'mb1','mb2','mb3','mb4','mb5',
            'tn1','tn2','tn3',
            'ss1','ss2','ss3'],
        ['jdf1','jdf2','jdf3','jdf4', # Jdf to SoG
            'sji1', 'sji2', 'sog1','sog2',
            'sog3','sog4'],
        ['jdf1','jdf2','jdf3','jdf4', # JdF to HC
            'ai1', 'ai2', 'ai3',
            'hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
        ['jdf1','jdf2','jdf3','jdf4', # JdF to Whidbey
            'ai1', 'ai2', 'ai3', 'ai4',
            'wb1','wb2','wb3','wb4']
            ]:

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
        fn = indir + sect_name + '.p'
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
    elif cc == 1:
        lstr = 'JdF to Strait of Georgia'
    elif cc == 2:
        lstr = 'JdF to Hood Canal'
    elif cc == 3:
        lstr = 'JdF to Whidbey Basin'
        
    ThalMean[lstr] = (sect_list, qin, qout, qsin, qsout, sin, sout, dist)
    
    cc += 1

# save results for plotting
pickle.dump(ThalMean, open(indir + 'ThalMean.p', 'wb'))


