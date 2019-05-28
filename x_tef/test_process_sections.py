# -*- coding: utf-8 -*-
"""
Code to test the bulk_calc code.
"""

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

import tef_fun_lorenz as tfl
from importlib import reload
reload(tfl)

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'

ex_list = ['cas6_v1_lo8_2017.01.01_2017.12.31', 'cas6_v2_lo8_2016.12.15_2017.12.31']
dd_offset_list = [0, 17*24]
dd_offset_dict = dict(zip(ex_list, dd_offset_list))

for ex in ex_list:
    fn = indir0 + ex + '/extractions/ai1.nc'
    dd_offset = dd_offset_dict[ex]

    print('\nWorking on ' + ex)
    
    # load fields
    ds = nc.Dataset(fn)
    q = ds['q'][dd_offset:, :, :]
    s = ds['salt'][dd_offset:, :, :]
    ot = ds['ocean_time'][dd_offset:]
    ds.close()
    
    print(Lfun.modtime_to_datetime(ot[0]))
    print(Lfun.modtime_to_datetime(ot[-1]))
    

    # TEF sort into salinity bins
    NT, NZ, NX = q.shape
    
    NT = 10
    
    
    # initialize intermediate results arrays for TEF quantities
    sedges = np.linspace(0, 36, 1001) # original was 1001 used 5001 for Willapa
    sbins = sedges[:-1] + np.diff(sedges)/2
    NS = len(sbins) # number of salinity bins

    # TEF variables
    tef_q = np.zeros((NT, NS))

    # other variables
    qnet = np.zeros(NT)

    for tt in range(NT):
        if np.mod(tt,100) == 0:
            print('  time %d out of %d' % (tt,NT))
            sys.stdout.flush()
        si = s[tt,:,:].squeeze()
        if isinstance(si, np.ma.MaskedArray):
            sf = si[si.mask==False].data.flatten()
        else:
            sf = si.flatten()
        qi = q[tt,:,:].squeeze()
        if isinstance(qi, np.ma.MaskedArray):
            qf = qi[si.mask==False].data.flatten()
        else:
            qf = qi.flatten()
        # sort into salinity bins
        inds = np.digitize(sf, sedges, right=True)
        counter = 0
        indsf = inds.copy().flatten()
        for ii in indsf:
            tef_q[tt, ii-1] += qf[counter]
            counter += 1
        
        # also keep track of volume transport
        qnet[tt] = qf.sum()
        
        print('\ntt = ' + str(tt))
        print(tef_q[tt,:].sum())
        print(qnet[tt].sum())
        
