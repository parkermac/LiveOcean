"""
Plot the mean of many TEF extractions on a thalweg section.

We organize by calling the saltier of the two layers "1" and the fresher "2".
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

# ************ NEW ********************************
# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir0)
indir0 = indir0 + item + '/'
indir = indir0 + 'bulk/'
outdir = indir0 + 'thalweg/'
Lfun.make_dir(outdir)
# **************************************************

# a structure to hold results for future use
ThalMean = dict()

# dict to define sections on each channel
channel_dict = {'JdF to South Sound':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3','ai4',
                'mb1','mb2','mb3','mb4','mb5',
                'tn1','tn2','tn3',
                'ss1','ss2','ss3'],
            'JdF to Strait of Georgia':['jdf1','jdf2','jdf3','jdf4',
                'sji1', 'sji2', 'sog1','sog2','sog3','sog4','sog5'],
            'JdF to Hood Canal':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3',
                'hc1','hc2','hc3','hc4','hc5','hc6','hc7','hc8'],
            'JdF to Whidbey Basin':['jdf1','jdf2','jdf3','jdf4',
                'ai1', 'ai2', 'ai3', 'ai4',
                'wb1','wb2','wb3','wb4','dp']}

for ch_str in channel_dict.keys():
    print('== ' + ch_str + ' ==')
    sect_list = channel_dict[ch_str]

    # initialize vectors to hold transports and salnities for all sections on a channel
    NS = len(sect_list)
    dd = np.nan * np.ones(NS)
    q1 = np.nan * np.ones(NS)
    q2 = np.nan * np.ones(NS)
    qs1 = np.nan * np.ones(NS)
    qs2 = np.nan * np.ones(NS)
    s1 = np.nan * np.ones(NS)
    s2 = np.nan * np.ones(NS)
    xs = np.nan * np.ones(NS)
    ys = np.nan * np.ones(NS)
    
    # then fill the vectors for that channel using the results of bulk_calc.py
    counter = 0
    for sect_name in sect_list:
        print(' ** ' + sect_name + ' **')
        
        # find x-y locations of the middle of each section
        x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
        lon = (x0+x1)/2
        lat = (y0+y1)/2    
        if counter == 0:
            lon0 = lon
            lat0 = lat
        xs[counter], ys[counter] = zfun.ll2xy(lon, lat, lon0, lat0)
        # also get section orientation
        if (x0==x1) and (y0!=y1):
            sdir = 'NS'
            a = [y0, y1]; a.sort()
            y0 = a[0]; y1 = a[1]
        elif (x0!=x1) and (y0==y1):
            sdir = 'EW'
            a = [x0, x1]; a.sort()
            x0 = a[0]; x1 = a[1]
            
        fn = indir + sect_name + '.p'
        bulk = pickle.load(open(fn, 'rb'))
        QQ = bulk['QQ']/1000 # convert to 1000 m3/s
        SS = bulk['SS']
        QQ1 = QQ.copy()
        QQ2 = QQ.copy()
        # initially we assume that the positive flux is the deeper (#1) layer
        QQ1[QQ<=0] = np.nan
        QQ2[QQ>0] = np.nan
        QQSS1 = QQ1 * SS
        QQSS2 = QQ2 * SS
        # here the nansum() is over all layers of a given sign from bulk_calc.py
        Q1 = np.nansum(QQ1,axis=1)
        Q2 = np.nansum(QQ2,axis=1)
        QS1 = np.nansum(QQSS1,axis=1)
        QS2 = np.nansum(QQSS2,axis=1)
        # here the nanmean() is averaging over all days
        # and the counter corresponds to which section this is on a channel
        q1[counter] = np.nanmean(Q1)
        q2[counter] = np.nanmean(Q2)
        qs1[counter] = np.nanmean(QS1)
        qs2[counter] = np.nanmean(QS2)
        s1[counter] = qs1[counter]/q1[counter]
        s2[counter] = qs2[counter]/q2[counter]
        # renumber when we got the direction wrong
        if s1[counter] < s2[counter]:
            print('   -- renumbering ---')
            # this cute trick with tuple unpacking reverses names 2,1 => 1,2
            q1[counter], q2[counter] = (q2[counter], q1[counter])
            qs1[counter], qs2[counter] = (qs2[counter], qs1[counter])
            s1[counter], s2[counter] = (s2[counter], s1[counter])
        else:
            pass
        counter += 1
        
    # create a distance vector (km)
    dx = np.diff(xs)
    dy = np.diff(ys)
    dd = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(NS)
    dist[1:] = np.cumsum(dd/1000)
    
    # pack results in a tuple as a dict entry
    ThalMean[ch_str] = (sect_list, q1, q2, qs1, qs2, s1, s2, dist)

# save results for plotting
pickle.dump(ThalMean, open(outdir + 'ThalMean.p', 'wb'))


