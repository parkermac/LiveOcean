"""
Process a TEF extraction.
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
indir = Ldir['LOo'] + 'extract/'

# choose the mooring extraction to process
print('\n%s\n' % '** Choose TEF extraction file to process **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.nc' in m) and ('tef_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if True:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
tef_file = m_dict[my_npt]
fn = indir + tef_file


# name output file
out_fn = indir + tef_file.replace('.nc','.p')
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
sdir = ds.sdir
ds0 = ds.date_string0
ds1 = ds.date_string1
ds.close()

# TEF calculation
qs = q*s
NT, NZ, NX = q.shape
# initialize intermediate results arrays for TEF quantities
sedges = np.linspace(0, 35, 1001) # original was 1001
sbins = sedges[:-1] + np.diff(sedges)/2
NS = len(sbins) # number of salinity bins
tef_q = np.zeros((NT, NS))
tef_qs = np.zeros((NT, NS))

# TEF variables
# which are also area integrals at ocean and river ends
tef_q = np.zeros((NT, NS))
tef_qs = np.zeros((NT, NS))

# other variables
qnet = np.zeros(NT)
fnet = np.zeros(NT)
g = 9.8
rho = 1025

for tt in range(NT): # NT
    if np.mod(tt,100) == 0:
        print('- time %d out of %d' % (tt,NT))
        sys.stdout.flush()
    si = s[tt,:,:].squeeze()
    sf = si[si.mask==False] # flattens the array
    qi = q[tt,:,:].squeeze()
    qf = qi[qi.mask==False]
    qsi = qs[tt,:,:].squeeze()
    qsf = qsi[qsi.mask==False]
    # sort into salinity bins
    inds = np.digitize(sf, sedges, right=True)
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


