"""
Plots mooring records.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os; import sys
sys.path.append(os.path.abspath('../alpha'))

import Lfun
import zfun


Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'moor/'

if False:
    # choose the mooring extraction to plot
    item = Lfun.choose_item(indir0)
    indir = indir0 + item + '/'
    infile = Lfun.choose_item(indir, tag='.nc')
else:
    mname = 'cas6_v3_lo8b_2019.05.01_2019.06.30'
    fname = 'Lahr1_hourly.nc'
fn = indir0 + mname + '/' + fname

ds = nc.Dataset(fn)

# time
ot_vec = ds['ocean_time'][:].data
mdate_vec = Lfun.modtime_to_mdate_vec(ot_vec)
mdt = mdates.num2date(mdate_vec) # list of datetimes of data

# space
zr =  ds['z_rho'][:]
Zr = zr.mean(axis=0)

# NOTE: velocities are packed (t,z)
# so we transpose when plotting using pcolormesh
u = ds['u'][:]
v = ds['v'][:]

if True:
    # Godin filter
    ulp = zfun.filt_godin_mat(u)
    vlp = zfun.filt_godin_mat(v)
    scl = .2
else:
    # Shorter time Hanning Filter
    ulp = zfun.filt_hanning_mat(u, n=6)
    vlp = zfun.filt_hanning_mat(v, n=6)
    scl = .1
    

up = u - ulp
vp = v - vlp

# rotate
Up = up[:,0]
Vp = vp[:,0]
th = 0.5 * np.arctan2(2*np.nanmean(Up*Vp),(np.nanvar(Up)-np.nanvar(Vp)))
cth = np.cos(th)
sth = np.sin(th)
urp = cth*up + sth*vp
vrp = cth*vp - sth*up

# plotting
plt.close('all')
fs=14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(18,10))

cmap = 'coolwarm'

dt0 = datetime(2019,5,26)
dt1 = datetime(2019,6,5)

ax = fig.add_subplot(211)
cs = ax.pcolormesh(mdt,Zr, urp.T, cmap=cmap, vmin = -scl, vmax=scl)
fig.colorbar(cs, ax=ax)
ax.set_xlim(dt0,dt1)
ax.grid(True)
#
ax.text(.05,.1,r'High-Passed Velocity [$m\ s^{-1}$] along $%d^{\circ}$ ($0^{\circ}=East$)' % (np.rad2deg(th)),
    transform=ax.transAxes)
#
ax.set_xticklabels([])
ax.set_ylabel('Z [m]')
ax.set_title(mname + ' ' + fname)

ax = fig.add_subplot(212)
cs = ax.pcolormesh(mdt,Zr, vrp.T, cmap=cmap, vmin = -scl, vmax=scl)
fig.colorbar(cs, ax=ax)
ax.set_xlim(dt0,dt1)
ax.grid(True)
#
ax.text(.05,.1,r'High-Passed Velocity [$m\ s^{-1}$] normal to $%d^{\circ}$' % (np.rad2deg(th)),
    transform=ax.transAxes)
#
ax.xaxis.set_major_locator(mdates.DayLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y.%m.%d'))
ax.xaxis.set_tick_params(labelrotation=25)
ax.set_xlabel('Date [UTC]')
ax.set_ylabel('Z [m]')

plt.show()
plt.rcdefaults()


