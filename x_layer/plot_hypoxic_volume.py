"""
Plots hypoxic volume time series.
"""


import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../plotting'))
import pfun

import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmocean
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'layer/'

# choose the layer extraction to plot
item = Lfun.choose_item(indir0)
indir = indir0 + item + '/'
fn = indir + 'hypoxic_volume_daily.nc'

# PLOTTING
#cmap = cmocean.cm.haline_r
cmap = cm.jet_r

plt.close('all')
fs=16
plt.rc('font', size=fs)
fig = plt.figure(figsize=(18,10))

ds = nc.Dataset(fn)

ot = ds['ocean_time'][:]
xp = ds['lon_psi'][:]; yp = ds['lat_psi'][:]
xr = ds['lon_rho'][:]; yr = ds['lat_rho'][:]
h = ds['h'][:]
h1 = np.ones_like(h)
mr = ds['mask_rho'][:]
DA = ds['DA'][:]
hyp_dz0 = ds['hyp_dz'][0,:,:]

# map
ax = fig.add_subplot(131)
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis([-130, -122, 42, 52])
ax.contour(xr,yr,h, [100, 200, 2000],
    colors=['darkorange','plum','darkorchid'], linewidths=2, linestyles='solid')
ax.text(.95,.16,'100 m',color='darkorange',weight='bold',transform=ax.transAxes,ha='right')
ax.text(.95,.13,'200 m',color='plum',weight='bold',transform=ax.transAxes,ha='right')
ax.text(.95,.1,'2000 m',color='darkorchid',weight='bold',transform=ax.transAxes,ha='right')
ax.set_title('Areas of Calculation')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks([-130, -128, -126, -124, -122])
ax.grid(True)

# time series
NT = len(ot)
days = (ot - ot[0])/86400.
mdays = Lfun.modtime_to_mdate_vec(ot)
mdt = mdates.num2date(mdays) # list of datetimes of data

dmax = 150
mask_dict = {}
mask_dict[49] = (yr >=49) & (yr < 50) & (xr < -125) & (h <= dmax) & (mr == 1)
mask_dict[48] = (yr >=48) & (yr < 49) & (xr < -124.6) & (h <= dmax) & (mr == 1)
mask_dict[47] = (yr >=47) & (yr < 48) & (xr < -123.5) & (h <= dmax) & (mr == 1)
mask_dict[46] = (yr >=46) & (yr < 47) & (xr < -123.5) & (h <= dmax) & (mr == 1)
mask_dict[45] = (yr >=45) & (yr < 46) & (xr < -123.5) & (h <= dmax) & (mr == 1)
mask_dict[44] = (yr >=44) & (yr < 45) & (xr < -123.5) & (h <= dmax) & (mr == 1)
mask_dict[43] = (yr >=43) & (yr < 44) & (xr < -124) & (h <= dmax) & (mr == 1)
NM = len(mask_dict)

lat_list = [49, 48, 47, 46, 45, 44, 43]

vol_dict = {}
for mm in lat_list:
    hm = np.ma.masked_where(~mask_dict[mm],h)
    h1m = np.ma.masked_where(~mask_dict[mm],h1)
    vol_dict[mm] = (hm*DA).sum()/1e9
    #print('lat = %d: area = %d [km3]' % (mm, vol_dict[mm]))
    ax.pcolormesh(xp,yp,mm*h1m[1:-1,1:-1], vmin=43, vmax=49, cmap=cmap)

hyp_vol = np.nan*np.ones((NT,NM))
hyp_rel_vol = np.nan*np.ones((NT,NM))
for ii in range(NT):
    if np.mod(ii,100) == 0:
        print('collecting ' + str(ii) + ' out of ' + str(NT))
        sys.stdout.flush()
    this_hdv = ds['hyp_dz'][ii,:,:] * DA
    mcount = 0
    for mm in lat_list:
        this_hdv_net = this_hdv[mask_dict[mm]].sum()/1e9
        hyp_vol[ii,mcount] = this_hdv_net
        hyp_rel_vol[ii,mcount] = this_hdv_net/vol_dict[mm]
        mcount += 1

ds.close()

hyp_vol_df = pd.DataFrame(hyp_vol,index=mdt, columns = lat_list)
hyp_rel_vol_df = pd.DataFrame(hyp_rel_vol,index=mdt, columns = lat_list)

ax1 = plt.subplot2grid((2,3), (0,1), colspan=2)
ax2 = plt.subplot2grid((2,3), (1,1), colspan=2)
for mm in lat_list:
    lc = cmap(int(255*(mm-43)/(49-43)))
    #print(lc)
    hyp_vol_df[mm].plot(ax=ax1, color=lc, legend=False, linewidth=2)
    hyp_rel_vol_df[mm].plot(ax=ax2, color=lc, legend=False, linewidth=2)
ax1.set_xticklabels([])
ax1.text(.03,.9, 'Hypoxic Volume [km3]', transform=ax1.transAxes)
ax2.text(.03,.9, 'Fractional Hypoxic Volume', transform=ax2.transAxes)

ax1.set_xlim(mdt[0], mdt[-1])
ax2.set_xlim(mdt[0], mdt[-1])

ax1.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])
ax2.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])

ax1.grid(True)
ax2.grid(True)

fig.tight_layout()
plt.show()
plt.rcdefaults()

