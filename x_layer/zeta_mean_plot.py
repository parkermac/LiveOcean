"""
Plots time-mean surface height field.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun

pth = os.path.abspath('../plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

Ldir = Lfun.Lstart()
fn = Ldir['LOo'] + 'layer/cas4_v2_lo6biom_2017.01.01_2017.12.31/zeta_hourly.nc'
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
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(xp,yp,zr[1:-1,1:-1], vmin=.13, vmax=.2)
fig.colorbar(cs)
pfun.dar(ax)
plt.show()
