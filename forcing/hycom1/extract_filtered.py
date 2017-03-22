#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:25:41 2017

@author: PM5

Extract filtered fields from LiveOcean "fh...p" files.

cd /Users/PM5
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname='cas1', tag='base')

import zfun

import pickle
import pandas as pd
from datetime import datetime, timedelta

# specify the output directory
out_dir = Ldir['data'] + 'hycom1/'

c = pickle.load(open(out_dir + 'coords_dict.p', 'rb'))
cc = c['91.2']

lon = cc['lon']
lat = cc['lat']
z = cc['z']

i0 = zfun.find_nearest_ind(lon, lon[0] + 1)
j0 = zfun.find_nearest_ind(lat, lat[0] + 1)
k0 = zfun.find_nearest_ind(z, -1500)

print('lon = %0.2f lat = %0.2f z = %0.2f' % (lon[i0], lat[j0], z[k0]))

in_dir0 = Ldir['LOo'] + Ldir['gtag'] + '/'

#%% generate a list of forecast days
dt0 = datetime(2013,1,1)
dt1 = datetime(2013,12,31)
dts_list = []
dt = dt0
while (dt1-dt).days > 0:
    dts = dt.strftime('%Y.%m.%d')
    dts_list.append(dts)
    dt = dt + timedelta(days=1)

fn_list = []
for dts in dts_list:
    fn_list.append(in_dir0 + 'f' + dts + '/ocn1/Data/fh' + dts + '.p')
 
#%% extract the data

DT = []
DTalt = []

vlist = ['ssh', 'u3d', 'v3d', 't3d', 's3d']

V = dict()
for v in vlist:
    V[v] = []

for fn in fn_list:
        
    b = pickle.load(open(fn, 'rb'))
    
    DT.append(b['dt']) # checking

    for v in vlist:
        if v == 'ssh':
            V[v].append(b[v][j0, i0])
        else:                
            V[v].append(b[v][k0, j0, i0])

#%% pack all into a dataframe
newind = pd.date_range(start=DT[0], end=DT[-1], freq='D') 
df = pd.DataFrame(columns = vlist, index=newind)

for v in vlist:
    df[v] = pd.Series(dict(zip(DT, V[v])), dtype=float)
    # for some reason I had to include the dtype or else
    # it did not parse the t3d and s3d as numbers
    
# and save for plotting
df.to_pickle(out_dir + 'extraction_fh_' + str(k0) + '.p')

