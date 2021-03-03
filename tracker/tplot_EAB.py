"""
Plot results of a particle tracking experiment.
Specific to a run from Elizabeth - trying to figure out how
deep particles from the slope get into JdF
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from time import time

Ldir = Lfun.Lstart()

# get Datasets
indir0 = '/Users/pm8/Documents/LiveOcean_output/tracks/'
indir = 'EAB_3d/'
rel = 'release_2017.09.24.nc'
fn = indir0 + indir + rel
fng = indir0 + indir + 'grid.nc'
dsr = nc.Dataset(fn)
dsg = nc.Dataset(fng)


# get a list of datetimes
ot_vec = dsr['ot'][:]
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp = dsg['lon_psi'][:]
latp = dsg['lat_psi'][:]
lonr = dsg['lon_rho'][:]
latr = dsg['lat_rho'][:]
hh = dsg['h'][:]
maskr = dsg['mask_rho'][:] # 1=water, 0=land

vdict = {}
vsdict = {}

tt0 = time()
for vn in ['lon', 'lat', 'z', 'cs', 'salt', 'temp']:
    vdict[vn] = dsr[vn][:]
print('time to load tracks = %0.2f sec' % (time()-tt0))

tt0 = time()
mask = (vdict['lon'][-1,:] > -124.5) & (vdict['lat'][-1,:] > 48) & (vdict['z'][0,:] < -500)
mask = zfun.fillit(mask)
mask = np.argwhere(mask).flatten()
print('time to make mask = %0.2f sec' % (time()-tt0))

tt0 = time()
for vn in vdict.keys():
    vsdict[vn] = vdict[vn][:,mask].copy()
print('time to subsample tracks = %0.2f sec' % (time()-tt0))

NT, NP = vdict['lon'].shape

NTS, NPS = vsdict['lon'].shape

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(16,10))

# MAP
# set domain limits
# automatically plot region of particles, with padding
pad = .1
aa = [vsdict['lon'].min() - pad, vsdict['lon'].max() + pad,
vsdict['lat'].min() - pad, vsdict['lat'].max() + pad]
ax = fig.add_subplot(121)
zm = -np.ma.masked_where(maskr==0, hh)
plt.pcolormesh(lonp, latp, zm[1:-1, 1:-1], vmin=-1000, vmax=0,
    cmap='terrain', alpha=.5)
plt.contour(lonr, latr, hh, [200, 300, 400, 500, 600, 700, 800, 900, 1000], colors='k')
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(indir.strip('/'))
ax.plot(vsdict['lon'], vsdict['lat'])

# time series
td = (ot_vec - ot_vec[0])/86400
tv_list = ['z', 'cs', 'salt', 'temp']
ntv = len(tv_list)
for ii in range(ntv):
    tv = tv_list[ii]
    NC = 2
    ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
    ax.plot(td, vsdict[tv])
    ax.text(.05, .05, tv, fontweight='bold', transform=ax.transAxes)
    if ii == ntv-1:
        ax.set_xlabel('Time (days)')

plt.show()

dsr.close()
dsg.close()

