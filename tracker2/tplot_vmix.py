"""
Plot results of a particle tracking experiment, specific to experiments about
vertical mixing of particles.

e.g. from:
python tracker.py -exp vmix -3d True -clb True -no_advection True
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

Ldir = Lfun.Lstart()

# get release Dataset
indir0 = Ldir['LOo'] + 'tracks2/'
indir = 'vmix_ndiv12_3d_nadv/'
rel = 'release_2019.07.04.nc'
dsr = nc4.Dataset(indir0 + indir + rel)

NT, NP = dsr['lon'].shape

# get a list of datetimes
ot_vec = dsr['ot'][:]
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]
t = (ot_vec - ot_vec[0])/3600

# Gather particle data
# packed [time, particle #]
lon = dsr['lon'][:]
lat = dsr['lat'][:]
z = dsr['z'][:]
h = dsr['h'][:]
salt = dsr['salt'][:]
temp = dsr['temp'][:]
cs = dsr['cs'][:]
zeta = dsr['zeta'][:]
dsr.close()

# rescale z to remove tides
ZZ = cs*h

# PLOTTING
plt.close('all')
fs = 14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(20,10))

# Histograms
for jj in [1,2,3]:
    
    zz = ZZ[:,1000*(jj-1):1000*jj - 1]
    ax = fig.add_subplot(1,3,jj)
    bins=np.linspace(zz[1,:].min(), 0, 50)
    for ii in range(1,NT-1, int(NT/10)):
        counts, obins = np.histogram(zz[ii,:], bins=bins)
        ax.plot(counts/NP, bins[:-1],'-o', label='Hour = %d' % (t[ii]))
    ax.set_xlim(0,)
    ax.set_xlabel('Fraction')
    ax.set_ylabel('Z [m]')
    ax.legend()


plt.show()
plt.rcdefaults()


