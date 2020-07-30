"""
Code to make a nice plot for intro figures of the first LiveOcean paper.

"""

# imports

#%% setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun, zrfun
import pfun

sys.path.append(os.path.abspath('../x_moor'))
import moor_lists

sys.path.append(os.path.abspath('../../ptools/tide_obs_mod'))
import obsfun as tide_ofn

from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import cmocean

# output location (assume it exists)
out_dir = Ldir['parent'] + 'LiveOcean_NOTES/Figures Slides Movies/'

# plotting choices
fs = 18
fs2 = 0.8*fs
fs3 = 0.6*fs2
vn = 'temp'; vmin=10; vmax=18
cmap = cmocean.cm.thermal
#cmap = 'RdYlBu_r'#'cubehelix'#'ocean'
aaf = [-125.3, -122.1, 46.8, 50.3] # focus domain
station = 'SAR003'

# get map fields
dstr = '2019.07.04'
fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f'+dstr+'/ocean_his_0020.nc'
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G = True)
lon = G['lon_psi'][:]
lat = G['lat_psi'][:]
v = ds[vn][0,-1,1:-1,1:-1]

# get river time series
r_df = pd.DataFrame(index=pd.date_range(start='1/1/2017', end='12/31/2019'),
    columns=['fraser','skagit'])
for year in [2017, 2018, 2019]:
    rfn = Ldir['LOo'] + 'river/cas6_v3_'+str(year)+'.01.01_'+str(year)+'.12.31.p'
    rdf = pd.read_pickle(rfn)
    # convert to 1000 m3/s
    rdf = rdf/1000
    rdt0 = datetime(year,1,1)
    rdt1 = datetime(year,12,31)
    # and save selected river(s) to the three year DataFrame
    r_df.loc[rdt0:rdt1,'fraser'] = rdf.loc[rdt0:rdt1,'fraser']
    r_df.loc[rdt0:rdt1,'skagit'] = rdf.loc[rdt0:rdt1,'skagit']

# get mooring record
mfn = (Ldir['LOo'] +
    'moor/cas6_v3_lo8b_2017.01.01_2019.12.31/' +
    station + '_daily.nc')
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

# get all mooring locations
moor_sta_dict, junk, junk, junk = moor_lists.get_sta_dict('shelf_moorings_mccabe')
moor_df = pd.DataFrame(index = moor_sta_dict.keys(), columns=['Longitude','Latitude'])
for sn in moor_sta_dict.keys():
    x, y = moor_sta_dict[sn]
    moor_df.loc[sn,'Longitude'] = x
    moor_df.loc[sn,'Latitude'] = y

# PLOTTING
plt.rc('font', size=fs)
plt.close('all')
fig = plt.figure(figsize=(13,13))

# -------------- FULL MAP -----------------------------------------------

ax = plt.subplot2grid((3,2), (0,0), rowspan=2)
cs = ax.pcolormesh(lon[6:-6,6:],lat[6:-6,6:],v[6:-6,6:], cmap=cmap, vmin=vmin, vmax=vmax)
nudge_alpha = .1
ax.pcolormesh(lon[:,:6],lat[:,:6],v[:,:6], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[:6,:],lat[:6,:],v[:6,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[-6:,:],lat[-6:,:],v[-6:,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
pfun.add_bathy_contours(ax, ds, txt=False)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(ds))
pfun.dar(ax)

# Inset colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left', borderpad=3) 
cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')

ax.text(.08, .53, r'SST $[^\circ C]$', transform=ax.transAxes)
ax.text(.08, .03, dstr, size=fs2, transform=ax.transAxes, style='italic')

ax.set_xticks([-130, -126, -122])
ax.set_yticks([42, 44, 46, 48, 50, 52])

# add ticks for grid spacing
x = lon[0,::10]
y = lat[::10,0]
for xx in x:
    ax.plot([xx,xx],[42,42.12],'-k',alpha=.5)
for yy in y:
    ax.plot([-122.18, -122],[yy,yy],'-k',alpha=.5)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

pfun.draw_box(ax, aaf, linestyle='-', color='k', alpha=.3, linewidth=3, inset=0)

ax.text(.93,.97,'(a) Full Model Domain', ha='right', va='top', weight='bold',
    transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='None', alpha=.6))

ah = ax.text(-125.456,49.4768,'Vancouver\nIsland', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)
ax.text(-123.072,46.7866,'Washington', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)
ax.text(-122.996,44.5788,'Oregon', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)

# -------------- FOCUS MAP -----------------------------------------------

ax = ax = plt.subplot2grid((3,2), (0,1), rowspan=2)
ax.pcolormesh(lon,lat,v, cmap=cmap, vmin=vmin, vmax=vmax)
pfun.add_coast(ax)
ax.axis(aaf)
pfun.dar(ax)

ax.set_xticks([-125, -124, -123])
ax.set_yticks([47, 48, 49, 50])

# add ticks for grid spacing
x = lon[0,::4]
y = lat[::4,0]
for xx in x:
    ax.plot([xx,xx],[aaf[2],aaf[2]+.06],'-k',alpha=.5)
for yy in y:
    hh = ax.plot([aaf[1]-.08, aaf[1]],[yy,yy],'-k',alpha=.5)

ax.set_xlabel('Longitude')

ax.text(.93,.97,'(b) Salish Sea', ha='right', va='top', size=fs,
    transform=ax.transAxes, weight='bold')

# add labels
ax.text(-122.682,49.335,'Fraser\nRiver',size=fs2,
    style='italic',ha='center',va='center',rotation=0)
ax.text(-123.7,49.2528,'Strait of Georgia',size=fs2,
    style='italic',ha='center',va='center',rotation=-30)
ax.text(-123.5,48.28,'Strait of Juan de Fuca',size=fs2,
    style='italic',ha='center',va='center',rotation=0,
    color='w')
ax.text(-123.3,47.6143,'Puget\nSound',size=fs2,
    style='italic',ha='center',va='center',rotation=+55)
ax.text(-122.3,48.48,'Skagit\nRiver',size=fs3,
    style='italic',ha='center',va='center',
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
ax.text(-123.173,48.44,'Haro\nStrait',size=fs3,
    style='italic',ha='center',va='center',
    color='w')
    
ds.close()


# -------------- RIVER TIME SERIES --------------------------------------

dt0 = datetime(2017,1,1,0)
dt1 = datetime(2020,1,1,0)

# river colors
rc1 = 'cornflowerblue'
rc2 = 'olive'

ax = fig.add_subplot(313)

r_df.plot(y='fraser', ax=ax, legend=False,
    xlim=(dt0,dt1), ylim=(0,12), grid=False, linestyle='-', lw=3, color=rc1)
r_df.plot(y='skagit', ax=ax, legend=False,
    xlim=(dt0,dt1), ylim=(0,12), grid=False, linestyle='-', lw=3, color=rc2)

ax.set_xticks([])
ax.set_xticks([], minor=True)
ax.vlines([datetime(2018,1,1),datetime(2019,1,1)],0,15, alpha=.5)

ax.text(.18,.5,'Fraser',color=rc1, transform=ax.transAxes, weight='bold')
ax.text(.12,.14,'Skagit',color=rc2, transform=ax.transAxes, weight='bold')

ax.set_xticks([datetime(2017,1,1),datetime(2017,7,1),datetime(2018,1,1),
    datetime(2018,7,1),datetime(2019,1,1),datetime(2019,7,1),datetime(2019,12,31)])
ax.set_xticklabels(['','2017','','2018','','2019',''], rotation=0,
    fontdict={'horizontalalignment':'center'})

ax.text(.95, .8, '(c) River Flow [$1000\ m^{3}s^{-1}$]',
    transform=ax.transAxes, weight='bold', ha='right',
    bbox=dict(facecolor='w', edgecolor='None'))


fig.tight_layout()
fig.savefig(out_dir + 'Fig_1_NEW.png')

plt.show()

plt.rcdefaults()




