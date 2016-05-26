"""
Plots mooring records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import numpy as np
from importlib import reload
import Lfun
reload(Lfun)
#import zfun  # plotting functions

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if '.nc' in m:
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
fn = m_dict[my_npt]

ds = nc.Dataset(indir + fn)

V = dict()
v1_list = ['ocean_time']
v2_list = ['zeta','sustr','svstr','swrad','lwrad',
    'shflux','latent','sensible']
v3_list_rho = ['u',
                 'v',
                 'temp',
                 'salt',
                 'NO3',
                 'phytoplankton',
                 'zooplankton',
                 'detritus',
                 'Ldetritus',
                 'oxygen',
                 'TIC',
                 'alkalinity',
                 'CaCO3',
                 'rho']
for vv in ds.variables:
    V[vv] = ds.variables[vv][:]

ds.close()

plt.close()

NR = 5; NC = 5
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])

cc = 0
nmid = round(V['salt'].shape[0]/2)
for vn in v3_list_rho:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(days, V[vn][-1,:], '-r')
    ax.plot(days, V[vn][nmid,:],'-g')
    ax.plot(days, V[vn][0,:], '-b')
    ax.set_xlim(0, days.max())
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_title(vn)
    cc += 1
for vn in v2_list:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(days, V[vn])
    ax.set_xlim(0, days.max())
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_title(vn)
    cc += 1
plt.show()
