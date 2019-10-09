"""
Code to make a more informative plot of a mooring extraction, focused
on physical properties.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
import seawater as sw

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='sj0', tag='v0')

sys.path.append(os.path.abspath('../plotting'))
import pfun


# choose the mooring extraction to plot
indir0 = Ldir['LOo'] + 'moor/'
if False:
    item = Lfun.choose_item(indir0)
    indir = indir0 + item + '/'
    infile = Lfun.choose_item(indir, tag='.nc')
else:
    indir = indir0 + 'sj0_v0_lo8nest_2019.09.01_2019.10.01/'
    infile = 'lost_tripod_hourly.nc'
fn = indir + infile

ds = nc.Dataset(fn)
lon = ds['lon_rho'][:]
lat = ds['lat_rho'][:]

dsg = nc.Dataset(Ldir['grid'] + 'grid.nc')

# get fields
ubar = ds['ubar'][:]
vbar = ds['vbar'][:]
u = ds['u'][:]
v = ds['v'][:]
salt = ds['salt'][:]
temp = ds['temp'][:]
pdens = sw.dens0(salt,temp)
z_w = ds['z_w'][:]
z_rho = ds['z_rho'][:]
ot = ds['ocean_time'][:]
t_md = Lfun.modtime_to_mdate_vec(ot)
t_dt = []
for t in ot:
    t_dt.append(Lfun.modtime_to_datetime(t))
    
def rot_vec(u,v,theta):
    ur = u*np.cos(theta) + v*np.sin(theta)
    vr = v*np.cos(theta) - u*np.sin(theta)
    return ur, vr
    
def de_mean(u, v):
    up = u - np.nanmean(u)
    vp = v - np.nanmean(v)
    return up, vp
    
# find angle of principal axes
up, vp = de_mean(ubar, vbar)
# we mask out high speeds because they introduced a funny direction
spd_lim = 0.5
spd = np.sqrt(up**2 + vp**2)
up[spd>spd_lim] = np.nan
vp[spd>spd_lim] = np.nan
up, vp = de_mean(up, vp)
theta = 0.5 * np.arctan2(2*np.nanmean(up*vp),(np.nanvar(up)-np.nanvar(vp)))

# and rotate
ubar_r, vbar_r = rot_vec(ubar, vbar, theta)
u_r, v_r = rot_vec(u,v,theta)

# plotting
plt.close('all')
fig = plt.figure(figsize=(18,10))

# map
ax = fig.add_subplot(3,4,4)
ax.plot(lon,lat,'*r')
pad = .1
ax.axis([lon-pad, lon+pad, lat-pad, lat+pad])
pfun.add_coast(ax)
pfun.dar(ax)

# depth averaged velocity scatterplot
V = 1
cmap = 'Spectral_r'
cmap = 'jet'
ax = fig.add_subplot(3,4,8)
ax.axhline(); ax.axvline()
ax.plot(ubar,vbar,'ob', alpha=.2, markersize=3)
ax.plot(ubar_r, vbar_r,'o', color='gold', alpha=.5, markersize=3)
ax.axis([-V,V,-V,V])
ax.grid(True)
ax.set_aspect('equal')
ax.axhline(); ax.axvline()

ax = plt.subplot2grid((2,4), (0,0), colspan=3)
T = (t_md)*np.ones((31,1))
tt = t_md - t_md[0]
cs = ax.pcolormesh(t_md, z_w.T, u_r.T - ubar_r, cmap=cmap, vmin=-V, vmax=V)
fig.colorbar(cs, ax=ax)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
ax.text(.05, .9, 'Along Channel (m/s Ebb Positive)', transform=ax.transAxes)
ax.set_ylim(top=5)

ax = plt.subplot2grid((2,4), (1,0), colspan=3)
T = (t_md)*np.ones((31,1))
tt = t_md - t_md[0]
cs = ax.pcolormesh(t_md, z_w.T, v_r.T -vbar_r, cmap=cmap, vmin=-V, vmax=V)
fig.colorbar(cs, ax=ax)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
ax.text(.05, .9, 'Cross Channel (m/s NE Positive)', transform=ax.transAxes)
ax.set_ylim(top=5)
#fig.autofmt_xdate() 3 rotates labels






plt.show()