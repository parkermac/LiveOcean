"""
Test of volume conservation for TEF.

RESULT:
cas4_v2_lo6biom_2017.01.01_2017.12.31

Qtef = -6221.7 (m3/s), Qr = 6260.8 (m3/s)
Error = 39.07 (m3/s) of 0.62%

which seems fine to me given that we did not account for dV/dt and the
transports were calculated from hourly snapshots.

"""

# imports
import pickle
import pandas as pd
import numpy as np

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart('cas4', 'v2')

import tef_fun
from importlib import reload
reload(tef_fun)

# select TEF input directory
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'

# form net flow through the sum of two sections at the bounds
# of the Salish Sea, negative means outflow
Qtef = 0
Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(indir + 'jdf1.p')
Qtef += qnet_lp.mean()
Qi, Si, Fi, qnet_lp, fnet_lp, td = tef_fun.tef_integrals(indir + 'sog4.p')
Qtef -= qnet_lp.mean()
# NOTE td goes 2.5 to 361.5 (360 entries, yearday)

# specify the rivers that flow into the volume bounded by these sections
salish_riv_list = [
    'clowhom',
    'squamish',
    'fraser',
    'tsolum',
    'oyster',
    'englishman',
    'cowichan',
    'nanaimo',
    'nooksack',
    'samish',
    'sanjuan',
    'hoko',
    'elwha',
    'dungeness',
    'dosewallips',
    'duckabush',
    'hamma',
    'skokomish',
    'deschutes',
    'nisqually',
    'puyallup',
    'green',
    'cedar',
    'snohomish',
    'stillaguamish',
    'skagit',
]

# gather the mean flow for each river, for a specific model/year
# Qr will be positive
fnr = Ldir['gtag'] + '_2017.01.01_2017.12.31.p'
fn = Ldir['LOo'] + 'river/' + fnr
df_riv = pd.read_pickle(fn)
Qr = 0
for rr in salish_riv_list:
    qr = df_riv[rr]
    # the times in the index are at hour zero of each day (366 entries)
    qrr = qr.iloc[2:-3].values # limit to times bounding the TEF td 
    qrmid = qrr[:-1] + np.diff(qrr)/2 # interpolate to be on TEF td times (360 values)
    qrm = qrmid.mean() # form the mean
    #print('%s %d' % (rr, int(qrm))) # interesting info
    Qr += qrm
    
# RESULTS
print('')
print('Qtef = %0.1f (m3/s), Qr = %0.1f (m3/s)' % (Qtef, Qr))
print('Error = %0.2f (m3/s) of %0.2f%%' % (Qtef+Qr, 100*(Qtef+Qr)/Qr))
    

    