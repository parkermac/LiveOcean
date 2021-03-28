"""
Code to make a nice plot of JdF Canyon bathymetry.

"""

# imports

#%% setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun
Ldir = Lfun.Lstart()
import pfun
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

# output location (assume it exists)
out_dir = Ldir['parent'] + 'LiveOcean_NOTES/Figures Slides Movies/'

sys.path.append(os.path.abspath('../x_moor'))
import moor_lists

# get all mooring locations
moor_sta_dict, junk, junk, junk = moor_lists.get_sta_dict('shelf_moorings_mccabe')
moor_df = pd.DataFrame(index = moor_sta_dict.keys(), columns=['Longitude','Latitude'])
for sn in moor_sta_dict.keys():
    x, y = moor_sta_dict[sn]
    moor_df.loc[sn,'Longitude'] = x
    moor_df.loc[sn,'Latitude'] = y



# get map fields
dstr = '2019.07.04'
fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f'+dstr+'/ocean_his_0020.nc'
G = zrfun.get_basic_info(fn, only_G = True)
lon = G['lon_psi']
lat = G['lat_psi']
z = -G['h']

zz = np.ma.masked_where(G['mask_rho']==False, z)

# plotting
plt.close('all')
fs=18
plt.rc('font', size=fs)

fig = plt.figure(figsize=(14,12))

ax = fig.add_subplot(121)
cs = ax.pcolormesh(lon, lat, zz[1:-1,1:-1], vmin=-3000, vmax = 100, cmap='terrain')
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis([-130, -122, 42, 52])
ax.set_xticks([-128, -126, -124])
ax.set_yticks([44, 46, 48, 50])
# Inset colorbar
cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right')
cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
cb.set_ticks([-2000, -200])
for t in cb.ax.get_xticklabels():
     t.set_fontsize(.8*fs)
aa1 = [-125.5, -124.5, 47.5, 49]
pfun.draw_box(ax, aa1, linestyle='-', color='r', alpha=1, linewidth=3, inset=0)
     
ax = fig.add_subplot(122)
cs = ax.pcolormesh(lon, lat, zz[1:-1,1:-1], vmin=-500, vmax = 0, cmap='terrain')
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(aa1)
ax.set_xticks([-125.5, -125, -124.5])
ax.set_yticks([47.5, 48, 48.5, 49])
# Inset colorbar
cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right')
cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
cb.set_ticks([-500, -200, 0])
for t in cb.ax.get_xticklabels():
     t.set_fontsize(.8*fs)
cs = ax.contour(G['lon_rho'], G['lat_rho'], zz, [-300, -200, -100],
    colors='k', linestyles='-', linewidths=1)
    
moor_df.plot(x='Longitude', y='Latitude',marker='o',ls='',mec='k',
    c='orange', ms=10,ax=ax,legend=False)



plt.show()

plt.rcdefaults()
