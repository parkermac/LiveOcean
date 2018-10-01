"""
Process TEF extractions.

Currently only set up to do salt.

"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

Ldir = Lfun.Lstart()

indir = Ldir['LOo'] + 'tef/'
if True:
    print('\nSelect an Extraction to plot:\n')
    List = os.listdir(indir)
    List.sort()
    NL = len(List)
    Ldict = dict(zip(range(NL), List))
    for ii in range(NL):
        print(str(ii) + ': ' + List[ii])
    my_ii = int(input('-- Input number: '))
    Litem = Ldict[my_ii]
else:
    Litem = 'cas4_v2_lo6biom_2017.01.01_2017.12.31'
print('\nProcessing ' + Litem + '\n')
Indir = indir + Litem + '/'

LList_raw = os.listdir(indir + Litem)
LList_raw.sort()

if False: # process all .nc files
    LList = [item for item in LList_raw if ('.nc' in item)]
else: # override
    LList = ['sog2.nc', 'sog3.nc']

for tef_file in LList:
    print(tef_file)
    fn = Indir + tef_file

    # name output file
    out_fn = Indir + tef_file.replace('.nc','.p')
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass

    # load fields
    ds = nc.Dataset(fn)
    q = ds['q'][:]
    s = ds['salt'][:]
    ot = ds['ocean_time'][:]
    zeta = ds['zeta'][:]
    gtagex = ds.gtagex
    ds0 = ds.date_string0
    ds1 = ds.date_string1
    ds.close()

    # TEF sort into salinity bins
    qs = q*s
    NT, NZ, NX = q.shape
    # initialize intermediate results arrays for TEF quantities
    sedges = np.linspace(0, 36, 1001) # original was 1001 used 5001 for Willapa
    sbins = sedges[:-1] + np.diff(sedges)/2
    NS = len(sbins) # number of salinity bins
    tef_q = np.zeros((NT, NS))
    tef_qs = np.zeros((NT, NS))

    # TEF variables
    tef_q = np.zeros((NT, NS))
    tef_qs = np.zeros((NT, NS))

    # other variables
    qnet = np.zeros(NT)
    fnet = np.zeros(NT)
    g = 9.8
    rho = 1025

    for tt in range(NT):
        if np.mod(tt,100) == 0:
            print('  time %d out of %d' % (tt,NT))
            sys.stdout.flush()
        si = s[tt,:,:].squeeze()
        sf = si[si.mask==False] # flattens the array
        qi = q[tt,:,:].squeeze()
        qf = qi[qi.mask==False]
        qsi = qs[tt,:,:].squeeze()
        qsf = qsi[qsi.mask==False]
        # sort into salinity bins
        inds = np.digitize(sf, sedges, right=True)
        inds[inds>NS] = NS # catching an addo Willaa error
        counter = 0
        for ii in inds:
            tef_q[tt, ii-1] += qf[counter]
            tef_qs[tt, ii-1] += qsf[counter]
            counter += 1
        
        # also keep track of volume transport
        qnet[tt] = qf.sum()
        
        # and tidal energy flux
        zi = zeta[tt,:].squeeze()
        ff = zi.reshape((1,NX)) * qi
        fnet[tt] = g * rho * ff.sum()

    # save results
    tef_dict = dict()
    tef_dict['tef_q'] = tef_q
    tef_dict['tef_qs'] = tef_qs
    tef_dict['sbins'] = sbins
    tef_dict['ot'] = ot
    tef_dict['qnet'] = qnet
    tef_dict['fnet'] = fnet
    pickle.dump(tef_dict, open(out_fn, 'wb'))


