"""
Plot results of a particle tracking experiment, specific to experiments about
why water is exchanged across the Admiralty Inlet sill.
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
indir = 'ai0_ndiv1_3d/'
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
u = dsr['u'][:]
v = dsr['v'][:]
w = dsr['w'][:]
salt = dsr['salt'][:]
temp = dsr['temp'][:]
cs = dsr['cs'][:]
zeta = dsr['zeta'][:]

z = cs*h

dsr.close()

# Choose the "winners"
mask = (lon[-1,:] > -122.75) & (lat[-1,:] < 48.1) # into Main Basin or HC
#mask = lat[-1,:] > 48.5 # into SoG
#mask = lon[-1,:] < -123.5 # out JdF
NPM = mask.sum()

# PLOTTING
plt.close('all')
fs = 16
plt.rc('font', size=fs)
fig = plt.figure(figsize=(18,10))

# Map
ax = fig.add_subplot(121)
ax.plot(lon[0,:], lat[0,:], '.k', alpha=.1)
ax.plot(lon[0,mask], lat[0,mask], '.c')
ax.plot(lon[-1,mask], lat[-1,mask], '.b')
pfun.dar(ax)
pfun.add_coast(ax)
aa = [-124, -122.3, 47.7, 49.1]
ax.axis(aa)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks([-124, -123.5, -123, -122.5])
ax.set_yticks([48, 48.5, 49])

# Histogram
ax = fig.add_subplot(222)
#zall = z.flatten()
z0 = z[0,mask]
z1 = z[-1,mask]
bins=np.linspace(-150,0,15 + 1)
#
counts, obins = np.histogram(z0, bins=bins)
ax.hist(bins[:-1], bins, weights=counts/NPM,
    orientation='horizontal', rwidth=.7, color='c')
#
counts, obins = np.histogram(z1, bins=bins)
ax.hist(bins[:-1], bins, weights=counts/NPM,
    orientation='horizontal', rwidth=.4, color='b')
ax.set_xlabel('Fraction')
ax.set_ylabel('Z [m]')

# Theta-S plot
ax = fig.add_subplot(224)
ax.plot(salt[0,mask], temp[0,mask],'.c', alpha=.7)
ax.plot(salt[-1,mask], temp[-1,mask],'+b', alpha=.4)
ax.set_xlabel('Salinity [g/kg]')
ax.set_ylabel('Potential Temp. [degC]')
# overlay potential density contours
ss = np.linspace(28,34,100)
tt = np.linspace(4,18,100)
SS, TT = np.meshgrid(ss,tt)
RR = sw.dens0(SS,TT)
CS = ax.contour(SS, TT, RR-1000)
ax.clabel(CS, fmt='%d')

plt.show()
plt.rcdefaults()


