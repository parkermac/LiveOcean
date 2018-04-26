"""
Code to plot a single day of hycom fields processed
by ocn1/make_forcing_main.py
"""

import os
import sys
import pickle
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

import pfun

pth = os.path.abspath('../forcing/ocn1')
if pth not in sys.path:
    sys.path.append(pth)
import Ofun

pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cas3', tag='v1')
import zfun
import zrfun

date_string = '2017.01.01'
frc = 'ocn1'

Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
Ldir['LOogf'] = (Ldir['LOog'] + 'f' + date_string + '/')
Ldir['LOogf_f'] = (Ldir['LOogf'] + frc + '/')
Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')

# get the ROMS grid
fng = Ldir['grid'] + 'grid.nc'
dsg = nc.Dataset(fng)
xr = dsg['lon_rho'][:]
yr = dsg['lat_rho'][:]
mr = dsg['mask_rho'][:]

# get the ini field on the ROMS grid
fn = Ldir['LOogf_f'] + 'ocean_ini.nc'
ds = nc.Dataset(fn)
tr = ds['temp'][0, -1, :, :].squeeze()
# mask on land
tr[mr==0] = np.nan

# get the processed hycom data that went into it
in_dir = Ldir['LOogf_fd']
# coordinates
xyz = pickle.load(open(in_dir + 'coord_dict.p', 'rb'))
x = xyz['lon']
y = xyz['lat']
# fields
fh = pickle.load(open(in_dir + 'fh' + date_string +'.p', 'rb'))
xfh = pickle.load(open(in_dir + 'xfh' + date_string +'.p', 'rb'))
t = fh['t3d'][-1,:,:]
xt = xfh['t3d'][-1,:,:]

# try adding a point and then extrapolating
i0, i1, ifr = zfun.get_interpolant(np.array([-122.75]), x)
j0, j1, jfr = zfun.get_interpolant(np.array([47.4]), x)
tt = t.copy()
tt[j0,i0] = 30 # doing this sets the mask to False automatically
X, Y = np.meshgrid(x,y)
TT = Ofun.extrap_nearest_to_masked(X, Y, tt)

# plotting
plt.close('all')
fig, axes = plt.subplots(1,4, squeeze=False, figsize=(22,8))

v0 = 7
v1 = 14

ax = axes[0,0]
cs = ax.pcolormesh(x,y,t, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])

ax = axes[0,1]
cs = ax.pcolormesh(x,y,xt, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])

ax = axes[0,2]
cs = ax.pcolormesh(x,y,TT, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])

ax = axes[0,3]
cs = ax.pcolormesh(xr,yr,tr, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])


plt.show()


