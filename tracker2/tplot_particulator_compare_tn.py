"""
Plot results of a particle tracking experiment and compare results to
the identical particulator experiment.
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import seawater as sw

from scipy import io

Ldir = Lfun.Lstart()

# get release Dataset
indir0 = Ldir['LOo'] + 'tracks2/'

indir = 'tn_ndiv12_surf/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)
# Gather particle data
# packed [time, particle #]
lon = dsr['lon'][:]
lat = dsr['lat'][:]
dsr.close()


# get the particulator output
pp = io.loadmat(Ldir['parent'] + 'particulator_output/test_tn.mat')
lonp = pp['x']
latp = pp['y']

# PLOTTING
plt.close('all')
fs = 14
lw = .3
alpha = .4
ms = 3
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,10))

pad = .1
aa = [lon.min()-pad, lon.max()+pad, lat.min()-pad, lat.max()+pad]

# Maps

ax = fig.add_subplot(121)
ax.plot(lon[0,:], lat[0,:],'.b', alpha=alpha, ms=ms)
ax.plot(lon, lat,'-', lw=lw, alpha=alpha)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('tracker2')

ax = fig.add_subplot(122)
ax.plot(lonp[0,:], latp[0,:],'.r', alpha=alpha, ms=ms)
ax.plot(lonp, latp,'-', lw=lw, alpha = alpha)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('particulator')


plt.show()
plt.rcdefaults()


