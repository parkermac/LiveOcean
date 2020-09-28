"""
Plots layer records for a "sandlance" extraction.

Designed to be standalone code for Matt Baker FHL Reesearch Apprenticeship.

2020_09 Parker MacCready
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

# USER edit path to reflect where the NetCDF file is
indir0 = '/Users/pm8/Documents/LiveOcean_output/layer/'
indir =indir0 + 'sj0_v0_lo8nest_2019.08.01_2019.08.10/'
fn = indir + 'sandlance_hourly.nc'

ds = nc.Dataset(fn)

# define convenience functions
def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))
def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx
    
# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(20,8))

x = ds['lon_rho'][:]
y = ds['lat_rho'][:]
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]

# v_list can contain any of the extracted fields:
# salt, temp, u, v, bustr, bvstr, ubar, vbar,
v_list = ['salt','temp','ubar', 'vbar']

F = dict() # field
S = dict() # series
ii = find_nearest_ind(x[0,:], -122.964)
jj = find_nearest_ind(y[:,0], 48.5)
for vn in v_list:
    F[vn] = ds[vn][0,:,:].squeeze()
    S[vn] = ds[vn][:,jj, ii].squeeze()

days = (ot - ot[0])/86400.

NC = len(v_list)
count = 1
for vn in v_list:
    ax = fig.add_subplot(2,NC,count)
    cs = ax.pcolormesh(xp, yp, F[vn][1:-1, 1:-1], cmap='Spectral')
    fig.colorbar(cs, ax=ax)
    ax.set_title(vn)
    dar(ax)
    ax.plot(x[jj,ii], y[jj,ii], '*r')
    count += 1

for vn in v_list:
    ax = fig.add_subplot(2,NC,count)
    ax.plot(days, S[vn])
    ax.text(.03, .9, vn, transform=ax.transAxes)
    try:
        ax.text(.03, .85, ds[vn].units, transform=ax.transAxes)
    except AttributeError:
        pass
    ax.text(.03, .8, ds[vn].long_name, transform=ax.transAxes)
    count += 1

plt.suptitle(fn)

plt.show()