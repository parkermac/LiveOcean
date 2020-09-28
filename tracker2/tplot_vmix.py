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

if True:
    # Choose an experiment to plot from.
    indir0 = Ldir['LOo'] + 'tracks2/'
    indir_list_raw = os.listdir(indir0)
    indir_list = []
    for d in indir_list_raw:
        if os.path.isdir(indir0 + d):
            indir_list.append(d)
    indir_list.sort()
    Npt = len(indir_list)#
    print('\n%s\n' % '** Choose Experiment to plot **')
    for npt in range(Npt):
        print(str(npt) + ': ' + indir_list[npt])
    my_npt = input('-- Experiment number (return = 0) --')
    if len(my_npt)==0:
        my_npt = 0
    indir = indir_list[int(my_npt)] + '/'

    # Choose a release from this experiment.
    rel_list = [rel for rel in os.listdir(indir0 + indir) if 'release' in rel]
    rel_list.sort()
    Nrl = len(rel_list)
    print('\n%s\n' % '** Choose Release file to plot **')
    for nrl in range(Nrl):
        print(str(nrl) + ': ' + rel_list[nrl])
    my_nrl = input('-- Release number (return = 0) -- ')
    if len(my_nrl)==0:
        my_nrl = 0
    rel = rel_list[int(my_nrl)]
else:
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
title_list = ['Slope', 'Juan de Fuca', 'Whidbey Basin']
for jj in [1,2,3]:
    
    NN = int(lon.shape[1]/3)
    zz = ZZ[:,NN*(jj-1):NN*jj - 1]
    ax = fig.add_subplot(1,3,jj)
    bins=np.linspace(zz[1,:].min(), 0, 100)
    for ii in range(0,NT-1, int(NT/10)):
        counts, obins = np.histogram(zz[ii,:], bins=bins)
        ax.plot(counts/NP, bins[:-1],'-o', label='Hour = %d' % (t[ii]))
    ax.set_xlim(0,0.01)
    ax.set_xlabel('Fraction')
    ax.set_ylabel('Z [m]')
    if jj==1:
        ax.legend()
    ax.set_title(title_list[jj-1])

plt.show()
plt.rcdefaults()


