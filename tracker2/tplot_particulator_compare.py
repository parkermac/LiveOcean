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

indir = 'jdf0_ndiv12_surf/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)
# Gather particle data
# packed [time, particle #]
lon12 = dsr['lon'][:]
lat12 = dsr['lat'][:]
dsr.close()

indir = 'jdf0_ndiv1_surf/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)
# Gather particle data
# packed [time, particle #]
lon1 = dsr['lon'][:]
lat1 = dsr['lat'][:]
dsr.close()


# get the particulator output
pp = io.loadmat(Ldir['parent'] + 'particulator_output/test3.mat')
lonp = pp['x']
latp = pp['y']

# PLOTTING
plt.close('all')
fs = 14
lw = .5
plt.rc('font', size=fs)
fig = plt.figure(figsize=(22,10))

pad = .1
aa = [lon1.min()-pad, lon1.max()+pad, lat1.min()-pad, lat1.max()+pad]

# Maps

ax = fig.add_subplot(131)
ax.plot(lon1[0,:], lat1[0,:],'.b', label='tracker ndiv=1')
ax.plot(lonp[0,:], latp[0,:],'.r', label='particulator')
ax.legend(loc='lower left')
ax.plot(lon1, lat1,'-b', lw=lw)
ax.plot(lonp, latp,'-r', lw=lw)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)

ax = fig.add_subplot(132)
ax.plot(lon12[0,:], lat12[0,:],'.c', label='tracker ndiv=12')
ax.plot(lonp[0,:], latp[0,:],'.r', label='particulator')
ax.legend(loc='lower left')
ax.plot(lon12, lat12,'-c', lw=lw)
ax.plot(lonp, latp,'-r', lw=lw)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)


ax = fig.add_subplot(133)
ax.plot(lon1[0,:], lat1[0,:],'.b', label='tracker ndiv=1')
ax.plot(lon12[0,:], lat12[0,:],'.c', label='tracker ndiv=12')
ax.legend(loc='lower left')
ax.plot(lon1, lat1,'-b', lw=lw)
ax.plot(lon12, lat12,'-c', lw=lw)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)


plt.show()
plt.rcdefaults()


