#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 09:08:41 2017

@author: PM5

Extract raw fields from the hycom1 archive.
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

import zfun

import pickle
import pandas as pd

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

a = os.listdir(out_dir)

aa = [i for i in a if i[0]=='h']

A = aa#[:300] # use [] to limit range for testing

DT = []
DTalt = []

vlist = ['ssh', 'u3d', 'v3d', 't3d', 's3d']

V = dict()
for v in vlist:
    V[v] = []

for fn in A:
    
#    dts = fn.strip('h').strip('.p')
#    dt = datetime.strptime(dts, '%Y.%m.%d')
#    DT.append(dt)
    
    b = pickle.load(open(out_dir + fn, 'rb'))
    
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
df.to_pickle(out_dir + 'extraction_raw_' + str(k0) + '.p')

