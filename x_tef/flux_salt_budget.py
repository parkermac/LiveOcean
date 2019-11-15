"""
Plot a basic salt budget for the segments.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import netCDF4 as nc
import argparse
from datetime import datetime, timedelta

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun

import tef_fun
import flux_fun

from time import time

# Get Ldir
Ldir = Lfun.Lstart('cas6', 'v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'

# select input/output location
indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
indir = indir0 + 'flux/'

# load DataFrames of transport, volume, and salinity
v_df = pd.read_pickle(indir + 'volumes.p')
s_df = pd.read_pickle(indir + 'daily_segment_salinity.p')

# seg list = list(v_df.index)
seg_list = flux_fun.seg_dict['Admiralty Inlet to South Sound']
seg_list.pop(0) # omit J4

vol = v_df.loc[seg_list,'volume m3'].sum()/1e9
print('volume = %0.1f km3' % (vol))

# form a time series of the net salt in the system
net_s_df = pd.DataFrame(index=s_df.index, columns=['net_salt'])
for dt in s_df.index:
    net_s_df.loc[dt,'net_salt'] = (s_df.loc[dt,seg_list]*v_df.loc[seg_list,'volume m3']).sum()
