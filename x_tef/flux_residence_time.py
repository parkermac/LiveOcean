"""
Designed to calculate residence time from all of the "IC" experiments.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'
import zfun

import tef_fun
import flux_fun

# Input directory
indir0 = Ldir['LOo'] + 'tef/'
indir = indir0 + 'flux_engine/' + Ldir['gtagex'] + '/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
V = flux_fun.get_V(v_df)

ic_list = [item for item in os.listdir(indir) if item[:3]=='IC_']

for infile in ic_list:
    source = infile.split('_')[0] + '_' + infile.split('_')[1]
    
    seg2_list = flux_fun.ic_seg2_dict[source]
        
    aa = pd.read_pickle(indir + infile)
    
    this_aa = aa.loc[:,seg2_list]
    this_V = V[seg2_list]
    net_V = this_V.sum()

    this_net_aa = this_aa.copy()

    for sn in this_V.index:
        VV = this_V[sn]
        this_net_aa.loc[:,sn] = this_net_aa.loc[:,sn] * VV
    
    # make a series of mean concentration in the volume
    mean_c = pd.Series(0,index=this_aa.index)
    mean_c = this_net_aa.sum(axis=1) / net_V

    # find e-folding time
    td = mean_c.index.values
    mc = mean_c.values

    ind_ef = np.argwhere(mc < 1/np.e)[0]
    
    print(' %10s: Tres = %0.2f days' % (infile.replace('.p',''), td[ind_ef]))


