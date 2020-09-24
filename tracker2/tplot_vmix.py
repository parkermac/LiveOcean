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
zz = cs*h

# PLOTTING
plt.close('all')
fs = 16
plt.rc('font', size=fs)
fig = plt.figure(figsize=(18,10))

# Histograms
ax = fig.add_subplot(111)
z0 = zz[1,:]
z0 = zz[int(NT/4),:]
z0 = zz[int(NT/2),:]
z0 = zz[int(3*NT/4),:]
z1 = zz[-1,:]

bins=np.linspace(zz[1,:].min(), 0, 200)
#
for ii in [1, int(NT/4), int(NT/2), int(3*NT/4), NT-1]:
    counts, obins = np.histogram(zz[ii,:], bins=bins)
    ax.hist(bins[:-1], bins, weights=counts/NP,
        orientation='horizontal', rwidth=.7, alpha = .3)
        
ax.set_xlabel('Fraction')
ax.set_ylabel('Z [m]')


plt.show()
plt.rcdefaults()


