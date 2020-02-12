"""
Code to make a nice plot for Figure 1 of the first LiveOcean paper.

"""

# imports

#%% setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun, zrfun
import pfun

from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
import matplotlib.pyplot as plt

# plotting choices
lwt = .5 # line thickness for grid resolution ticks
fs = 13
#vn = 'salt'; vmin=20; vmax=33
vn = 'temp'; vmin=8; vmax=20
cmap = 'RdYlBu_r'#'cubehelix'#'ocean'
aaf = [-125.3, -122.1, 46.8, 50.3] # focus domain

# get map fields
fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f2019.07.04/ocean_his_0020.nc'
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G = True)
lon = G['lon_psi'][:]
lat = G['lat_psi'][:]
v = ds[vn][0,-1,1:-1,1:-1]

# get river time series
r_df = pd.DataFrame(index=pd.date_range(start='1/1/2017', end='12/31/2019'), columns=['fraser'])
for year in [2017, 2018, 2019]:
    rfn = Ldir['LOo'] + 'river/cas6_v3_'+str(year)+'.01.01_'+str(year)+'.12.31.p'
    rdf = pd.read_pickle(rfn)
    # convert to 1000 m3/s
    rdf = rdf/1000
    rdt0 = datetime(year,1,1)
    rdt1 = datetime(year,12,31)
    # and save selected river(s) to the three year DataFrame
    r_df.loc[rdt0:rdt1,'fraser'] = rdf.loc[rdt0:rdt1,'fraser']

# get mooring record
mfn = (Ldir['LOo'] +
    'moor/cas6_v3_lo8b_2017.01.01_2019.12.31/' +
    'HCB010_hourly.nc')
m_ds = nc.Dataset(mfn)
mlon = m_ds['lon_rho'][:]
mlat = m_ds['lat_rho'][:]
ot_vec = m_ds['ocean_time'][:]
dt_list = []
for ot in ot_vec:
    dt_list.append(Lfun.modtime_to_datetime(ot))
m_df = pd.DataFrame(index=dt_list, columns=['sustr', 'svstr','salt_top','salt_bot'])
m_df.loc[:,'sustr'] = m_ds['sustr'][:]
m_df.loc[:,'svstr'] = m_ds['svstr'][:]
m_df.loc[:,'salt_top'] = m_ds['salt'][:,-1] # top salinity
m_df.loc[:,'salt_bot'] = m_ds['salt'][:,0] # bottom salinity
# remove 1/1/2020
m_df = m_df.iloc[:-1,:]
# and make daily averages
md_df = m_df.resample('1D').mean()


# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(14,7))

# full map
ax = fig.add_subplot(131)
ax.pcolormesh(lon[6:-6,6:],lat[6:-6,6:],v[6:-6,6:], cmap=cmap, vmin=vmin, vmax=vmax)
nudge_alpha = .1
ax.pcolormesh(lon[:,:6],lat[:,:6],v[:,:6], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[:6,:],lat[:6,:],v[:6,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[-6:,:],lat[-6:,:],v[-6:,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
pfun.add_bathy_contours(ax, ds, txt=False)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(ds))
pfun.dar(ax)
ax.set_xticks([-130, -126, -122])
ax.set_yticks([42, 44, 46, 48, 50, 52])
ax.tick_params(labelsize=fs) # tick labels
# add ticks for grid spacing
x = lon[0,::10]
y = lat[::10,0]
for xx in x:
    ax.plot([xx,xx],[42,42.12],'-k', lw=lwt)
for yy in y:
    ax.plot([-122.18, -122],[yy,yy],'-k', lw=lwt)
ax.set_xlabel('Longitude (deg)', size=fs)
ax.set_ylabel('Latitude (deg)', size=fs)
pfun.draw_box(ax, aaf, linestyle='-', color='k', alpha=.3, linewidth=3, inset=0)
ax.set_title('(a) Full Model Domain', size=fs)
ah = ax.text(-126.27,48.65,'Vancouver\nIsland',size=fs-2,style='italic',rotation=-45)
ax.text(-124,46.1,'Washington',size=fs-2,style='italic',rotation=-45)
ax.text(-123.7,44.2,'Oregon',size=fs-2,style='italic',rotation=-45)
    
# focus map
ax = fig.add_subplot(132)
ax.pcolormesh(lon,lat,v, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.add_coast(ax)
ax.axis(aaf)
pfun.dar(ax)
ax.set_xticks([-125, -124, -123])
ax.set_yticks([47, 48, 49, 50])
ax.tick_params(labelsize=fs) # tick labels
# add ticks for grid spacing
x = lon[0,::4]
y = lat[::4,0]
for xx in x:
    ax.plot([xx,xx],[aaf[2],aaf[2]+.06],'-k', lw=lwt)
for yy in y:
    hh = ax.plot([aaf[1]-.08, aaf[1]],[yy,yy],'-k', lw=lwt)
ax.set_xlabel('Longitude (deg)', size=fs)
# add mooring location
ax.plot(mlon, mlat, 'sk', ms=4)
ax.set_title('(b) Salish Sea', size=fs)

dt0 = datetime(2017,1,1,0)
dt1 = datetime(2020,1,1,0)

# river time series
ax = fig.add_subplot(333)
hh = r_df.plot(y='fraser', label='Fraser R. (1000 m3/s)', ax=ax, legend=True,
    xlim=(dt0,dt1), ylim=(0,15), grid=True, style='-g')
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.text(.05, .8, '(c)', size=fs, transform=ax.transAxes)

# wind time series
ax = fig.add_subplot(336)
md_df.plot(y='svstr', label='N-S Windstress (Pa)', ax=ax, legend=True, xlim=(dt0,dt1), grid=True, style='-c')
ax.set_xticklabels([])
ax.set_xticklabels([], minor=True)
ax.text(.05, .1, '(d)', size=fs, transform=ax.transAxes)

# salinity time series
ax = fig.add_subplot(339)
md_df.plot(y=['salt_bot', 'salt_top'], label=['Bottom Salinity', 'Surface Salinity'],
    ax=ax, legend=True, xlim=(dt0,dt1), grid=True)
ax.text(.05, .1, '(e)', size=fs, transform=ax.transAxes)

#fig.tight_layout()
plt.show()

ds.close()


