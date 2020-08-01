"""
Plots layer records.

Customized for bottom DO for some events - for Sam and Emily 2020.07.31.

"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun

sys.path.append(os.path.abspath('../plotting'))
import pfun

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'layer/'

fn1 = indir0 + 'cas6_v3_lo8b_2017.07.17_2017.07.26/bottom_DO_daily.nc'
fn2 = indir0 + 'cas6_v3_lo8b_2018.06.07_2018.06.12/bottom_DO_daily.nc'

ds1 = nc.Dataset(fn1)
ds2 = nc.Dataset(fn2)

xp = ds1['lon_psi'][:]
yp = ds1['lat_psi'][:]
ot1 = ds1['ocean_time'][:]
ot2 = ds2['ocean_time'][:]

dt1a = Lfun.modtime_to_datetime(ot1[0])
dt1b = Lfun.modtime_to_datetime(ot1[-1])
dt2a = Lfun.modtime_to_datetime(ot2[0])
dt2b = Lfun.modtime_to_datetime(ot2[-1])

# PLOTTING
plt.close('all')
fs=14
plt.rc('font', size=fs)

fig = plt.figure(figsize=(16,10))

ax = fig.add_subplot(121)
f = ds1['oxygen'][:]
ff = f.mean(axis=0)
cs = ax.pcolormesh(xp,yp,ff[1:-1,1:-1],cmap=cmocean.cm.oxy,vmin=0,vmax=300)
pfun.dar(ax)
pfun.add_coast(ax)
pfun.add_bathy_contours(ax,ds1, txt=True)
ax.axis(pfun.get_aa(ds1))
fig.colorbar(cs, ax=ax)
ax.text(.03,.15,dt1a.strftime('%Y.%m.%d')+'\n      to\n'+dt1b.strftime('%Y.%m.%d'),
    weight='bold',c='w',transform=ax.transAxes)
ax.set_title(r'Bottom DO $[\mu M]$')

ax = fig.add_subplot(122)
f = ds2['oxygen'][:]
ff = f.mean(axis=0)
cs = ax.pcolormesh(xp,yp,ff[1:-1,1:-1],cmap=cmocean.cm.oxy,vmin=0,vmax=300)
pfun.dar(ax)
pfun.add_coast(ax)
pfun.add_bathy_contours(ax,ds1, txt=True)
ax.axis(pfun.get_aa(ds1))
fig.colorbar(cs, ax=ax)
ax.text(.03,.15,dt2a.strftime('%Y.%m.%d')+'\n      to\n'+dt2b.strftime('%Y.%m.%d'),
    weight='bold',c='w',transform=ax.transAxes)
ax.set_title(r'Bottom DO $[\mu M]$')


plt.show()

plt.rcdefaults()

