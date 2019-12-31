"""
Goes through all the sections for a given extraction and makes
long-time averages (like a whole year) and forces fluxes into
two layers.

The results are saved in a single DataFrame, for later use by code
like flux_get_A.py.

NEW: Now allowing for user-selected time periods.  Currently it saves the full year average
and two other "seasons".  See the dict "dtr".

"""

# setup
import numpy as np
import pickle
import pandas as pd
from datetime import datetime, timedelta

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun

import tef_fun
from importlib import reload
reload(tef_fun)

# user settings
testing = False

# desired time ranges
dtr = {}
dtr['full'] = (datetime(2017,1,1,12,0,0), datetime(2017,12,31,12,0,0))
dtr['spring'] = (datetime(2017,3,1,12,0,0), datetime(2017,6,1,12,0,0))
dtr['fall'] = (datetime(2017,9,1,12,0,0), datetime(2017,12,1,12,0,0))

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# choose input and organize output
Ldir = Lfun.Lstart()
indir00 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir00)
indir0 = indir00 + item + '/'
indir = indir0 + 'bulk/'

outdir = indir0 + 'flux/'
Lfun.make_dir(outdir)

if testing == True:
    sect_name_list = ['hc3','hc4']
else:
    sect_name_list = list(sect_df.index)

for season in dtr.keys():
    
    dt0 = dtr[season][0]
    dt1 = dtr[season][1]
    
    # initialize output DataFrame
    # naming convention part 1: q_, f_, s_ mean volume flux, salt flux, and salinity
    # naming convention part 2: _s and _f mean the saltier or fresher of the two layers
    # ** sign convention: positive Northward or Eastward **
    df = pd.DataFrame(index=sect_name_list, columns=['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'])
    
    get_time = True
    for sect_name in sect_name_list:
        print('== ' + sect_name + ' ==')
        # then fill the DataFrame for this section using the results of bulk_calc.py
        
        x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
        
        fn = indir + sect_name + '.p'
        bulk = pickle.load(open(fn, 'rb'))
    
        # find the index range that corresponds to the desired time range
        if get_time == True:
            ot_vec = bulk['ot']
            dt_list = []
            for ot in ot_vec:
                dt_list.append(Lfun.modtime_to_datetime(ot))
            get_time = False
        
        try:
            idt0 =dt_list.index(dt0)
        except ValueError:
            # assume we want to start at the begnning
            idt0 = 0
            if testing == True:
                print('** dt0 not found, using idt0 = %d' % (idt0))
        try:
            idt1 =dt_list.index(dt1)
        except ValueError:
            # assume we want to end at the last day
            idt1 = len(ot_vec)
            if testing == True:
                print('** dt1 not found, using idt1 = %d' % (idt1))
        
        QQ = landward * bulk['QQ'] # multiplying by "landward" satisfies the sign convention
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
        # put results into the naming convention
        # here the nanmean() is averaging over all days
        q_s = np.nanmean(Q1[idt0:idt1])
        q_f = np.nanmean(Q2[idt0:idt1])
        f_s = np.nanmean(QS1[idt0:idt1])
        f_f = np.nanmean(QS2[idt0:idt1])
        s_s = f_s/q_s
        s_f = f_f/q_f
        # renumber when we got the direction wrong
        if s_s < s_f:
            print('   -- resorting ---')
            # this cute trick with tuple unpacking reverses names _f, _s => _s, _f
            q_s, q_f = (q_f, q_s)
            f_s, f_f = (f_f, f_s)
            s_s, s_f = (s_f, s_s)
        else:
            pass
        # save results into the DataFrame
        df.loc[sect_name, 'q_s'] = q_s
        df.loc[sect_name, 'q_f'] = q_f
        df.loc[sect_name, 'f_s'] = f_s
        df.loc[sect_name, 'f_f'] = f_f
        df.loc[sect_name, 's_s'] = s_s
        df.loc[sect_name, 's_f'] = s_f
        df.loc[sect_name, 'lon'] = (x0 + x1)/2
        df.loc[sect_name, 'lat'] = (y0 + y1)/2
        
    if testing == False:
        # save results for plotting
        df.to_pickle(outdir + 'two_layer_' + season + '.p')
    else:
        print('========= ' + season + ' =============')
        print(df)


