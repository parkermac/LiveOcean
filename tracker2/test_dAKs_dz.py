"""
Code to test calculation of dAKs/dz for vmix.
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

fn = ('/Users/pm8/Documents/LiveOcean_output/moor/' +
    'cas6_v3_lo8b_2019.07.04_2019.07.04/vmixWB_hourly.nc')
ds = nc.Dataset(fn)


# plotting
plt.close('all')
fig = plt.figure(figsize=(12,12))

for ii in [24]:
    
    k = ds['AKs'][ii,:]
    zw = ds['z_w'][ii,:]
    zr = ds['z_rho'][ii,:]

    dz = np.diff(zw)
    dk = np.diff(k)
    dkdz = dk/dz

    ax = fig.add_subplot(121)
    ax.plot(k, zw, '-ob', lw=2)
    ax.set_title('AKs')
    ax.set_ylim(zw[0],zw[-1])

    ax = fig.add_subplot(122)
    ax.plot(dkdz, zr, '-or', lw=2)
    ax.set_title('dAKs/dz')
    ax.set_ylim(zw[0],zw[-1])

plt.show()