"""
Goes through all the sections for a given extraction and makes
long-time averages (like a whole year) and forces fluxes into
two layers.

The results are saved in a single DataFrame, for later use by code
like flux_get_A.py.

"""

# setup
import numpy as np
import pickle
import pandas as pd

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

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir0)
indir0 = indir0 + item + '/'
indir = indir0 + 'bulk/'
outdir = indir0 + 'flux/'
Lfun.make_dir(outdir)

# initialize output DataFrame
# naming convention part 1: q_, f_, s_ mean volume flux, salt flux, and salinity
# naming convention part 2: _s and _f mean the saltier or fresher of the two layers
# sign convention: positive Northward or Eastward
df = pd.DataFrame(index=sect_df.index, columns=['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'])

for sect_name in sect_df.index:
    print('== ' + sect_name + ' ==')
    # then fill the DataFrame for this section using the results of bulk_calc.py
        
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]    
        
    fn = indir + sect_name + '.p'
    bulk = pickle.load(open(fn, 'rb'))
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
    q_s = np.nanmean(Q1)
    q_f = np.nanmean(Q2)
    f_s = np.nanmean(QS1)
    f_f = np.nanmean(QS2)
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
        
# save results for plotting
df.to_pickle(outdir + 'two_layer.p')


