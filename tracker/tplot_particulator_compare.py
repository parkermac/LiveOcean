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
indir0 = Ldir['LOo'] + 'tracks/'

indir = 'eddy0_ndiv12_3d/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)
# Gather particle data
# packed [time, particle #]
lon = dsr['lon'][:]
lat = dsr['lat'][:]
h = dsr['h'][:]
zeta = dsr['zeta'][:]
cs = dsr['cs'][:]
dsr.close()

# rescale z to remove tides
#z = cs*(h+zeta)
z = cs*h

# get the particulator output
pp = io.loadmat(Ldir['parent'] + 'particulator_output/eddy0.mat')
lonp = pp['x']
latp = pp['y']
#zp = pp['sigma']*(pp['H']+pp['zeta'])
zp = pp['sigma']*pp['H']


# PLOTTING
plt.close('all')
fs = 14
lw = .3
alpha = .4
ms = 3
plt.rc('font', size=fs)
fig = plt.figure(figsize=(16,12))

pad = .1
aa = [lon.min()-pad, lon.max()+pad, lat.min()-pad, lat.max()+pad]

# Maps

ax = fig.add_subplot(221)
ax.plot(lon[0,:], lat[0,:],'.b', alpha=alpha, ms=ms)
ax.plot(lon, lat,'-', lw=lw, alpha=alpha)
ax.plot(lon[-1,:], lat[-1,:],'.b', ms=ms)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('tracker')

ax = fig.add_subplot(222)
ax.plot(lonp[0,:], latp[0,:],'.r', alpha=alpha, ms=ms)
ax.plot(lonp, latp,'-', lw=lw, alpha = alpha)
ax.plot(lonp[-1,:], latp[-1,:],'.r', ms=ms)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('particulator')

ax = fig.add_subplot(212)

zmin = min(z.min(), zp.min())

NT, NP = z.shape
NTp, NPp = zp.shape
bins=np.linspace(zmin, 0, 15)

for ii in range(1,10):
    counts, obins = np.histogram(z[-ii,:], bins=bins)
    ax.plot(counts/NP, bins[:-1],'-ob', alpha=alpha)
        
    counts, obins = np.histogram(zp[-ii,:], bins=bins)
    ax.plot(counts/NPp, bins[:-1],'-or', alpha=alpha)
        
ax.set_xlim(0,)
ax.set_ylim(-100,0)
ax.set_xlabel('Fraction')
ax.set_ylabel('Z [m]')


plt.show()
plt.rcdefaults()


