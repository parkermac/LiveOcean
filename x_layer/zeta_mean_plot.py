"""
Plots time-mean surface height field.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../plotting'))
import pfun

Ldir = Lfun.Lstart()
fn = Ldir['LOo'] + 'layer/cas6_v3_lo8b_2019.06.01_2019.08.31/surface_hourly.nc'
ds = nc.Dataset(fn)

# gather fields
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]
mr = ds['mask_rho'][:]
NY, NX = mr.shape
NT = len(ot)

# form mean zeta
zr = np.zeros_like(mr)

for tt in range(NT):
    if np.mod(tt,100)==0:
        print('Step %d out of %d' % (tt,NT))
    zr = zr + ds['zeta'][tt,:,:].squeeze()/NT

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(14,10))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(xp,yp,zr[1:-1,1:-1], vmin=-.1, vmax=.1, cmap='nipy_spectral')
ax.axis([-126, -122, 47, 50])
fig.colorbar(cs)
pfun.dar(ax)
plt.show()
