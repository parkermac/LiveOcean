#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:50:41 2017

@author: PM5

Code to test the routines that extrapolate HyCOM files.

"""

import os
import sys
fpth = os.path.abspath('../../alpha')
if fpth not in sys.path:
    sys.path.append(fpth)
import Lfun
Ldir = Lfun.Lstart('cas1', 'base')
Ldir['date_string'] = '2013.01.01'
Ldir['LOogf_fd'] = (Ldir['LOo'] + Ldir['gtag'] + '/f' + Ldir['date_string']
        + '/ocn1/Data/')

import matplotlib.pyplot as plt

import Ofun
from importlib import reload
reload(Ofun)

fh_dir = Ldir['LOogf_fd']

#%% make a list of files to work on

a = os.listdir(fh_dir)
aa = [item for item in a if item[:2]=='fh']
# but for this test we just work on one
fn = aa[1]
print(fn)

in_fn = fh_dir + fn

#%% do the extrapolated using function calls

lon, lat, z, L, M, N, X, Y = Ofun.get_coords(fh_dir)
V = Ofun.get_extrapolated(in_fn, L, M, N, X, Y, z)
     
#%% Plotting
plt.close('all')

# 2D fields
fig, axes = plt.subplots(1, 3, figsize=(13,8))
vn2_list = ['ssh', 'ubar', 'vbar']
for ii in range(3):
    vn = vn2_list[ii]
    ax = axes[ii]
    pc = ax.pcolormesh(lon, lat, V[vn])
    fig.colorbar(pc, ax=ax)
    ax.set_title(vn)

# 3D fields
vn3_list = ['theta', 's3d', 'u3d', 'v3d']
for vn in vn3_list:
    fig, axes = plt.subplots(3, 3, figsize=(13,8))    
    VV = V[vn]
    # maps at three depth levels    
    for ii in range(3):
        ax = axes[0, ii]
        nn = int(ii * (N-1)/2)
        pc = ax.pcolormesh(lon, lat, VV[nn, :, :])
        fig.colorbar(pc, ax=ax)
        if ii==0:
            ax.set_title(vn)
    # sections at three longitudes
    for ii in range(3):
        ax = axes[1, ii]
        ll = int(ii * (L-1)/2)
        pc = ax.pcolormesh(lat, z/1000, VV[:, :, ll])
        fig.colorbar(pc, ax=ax)
    # sections at three latitudes
    for ii in range(3):
        ax = axes[2, ii]
        mm = int(ii * (M-1)/2)
        pc = ax.pcolormesh(lon, z/1000, VV[:, mm, :])
        fig.colorbar(pc, ax=ax)
    
plt.show()
        