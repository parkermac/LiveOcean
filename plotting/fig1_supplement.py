"""
Code to make a nice plot for Supplement of the first LiveOcean paper, Figure 1.

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
fs3 = 0.75*fs2
vn = 'h'
vmin=0; vmax=3000
vmin2=0; vmax2=300

#cmap = cmocean.cm.thermal
cmap = cmocean.cm.deep #'RdYlBu_r'#'cubehelix'#'ocean'
#cmap = 'gist_earth_r'
#cmap='terrain_r'
#cmap='binary'
aaf = [-125.3, -122.1, 46.8, 50.3] # focus domain

# marker colors and styles
# [tide, ctd, mooring]
C = ['#1f77b4', '#ff7f0e', '#2ca02c']
M = ['d','*','o'] 

# get map fields
dstr = '2019.07.04'
fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f'+dstr+'/ocean_his_0020.nc'
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G = True)
lon = G['lon_psi'][:]
lat = G['lat_psi'][:]
if vn == 'h':
    v = ds[vn][1:-1,1:-1]
    mm = G['mask_rho'][1:-1,1:-1].data
    v[~mm] = np.nan
else:
    v = ds[vn][0,-1,1:-1,1:-1]

# get all mooring locations
moor_sta_dict, junk, junk, junk = moor_lists.get_sta_dict('shelf_moorings_mccabe')
moor_df = pd.DataFrame(index = moor_sta_dict.keys(), columns=['Longitude','Latitude'])
for sn in moor_sta_dict.keys():
    x, y = moor_sta_dict[sn]
    moor_df.loc[sn,'Longitude'] = x
    moor_df.loc[sn,'Latitude'] = y
        
# get all CTD station locations
sta_df = pd.read_pickle(Ldir['parent'] + 'ptools_data/ecology/' + 'sta_df.p')
sta_df_ca = pd.read_pickle(Ldir['parent'] + 'ptools_data/canada/' + 'sta_df.p')
sta_df = pd.concat((sta_df, sta_df_ca), sort=False)
ctd_df = sta_df.loc[:,['Longitude','Latitude']]

# get tide station locations
junk, junk, tide_sn_dict = tide_ofn.get_sn_dicts()
tide_df = pd.DataFrame(index = tide_sn_dict.keys(), columns=['Longitude','Latitude'])
for name in tide_sn_dict.keys():
    sn = tide_sn_dict[name] # station number
    mfn = Ldir['parent'] + 'ptools_output/tide/obs_data/m_' + str(sn) + '_' + str(2017) + '.csv'
    Mo = Lfun.csv_to_dict(mfn) # has name and lon,lat info for a station
    tide_df.loc[name,'Longitude'] = float(Mo['lon'])
    tide_df.loc[name,'Latitude'] = float(Mo['lat'])
    
def add_sta_loc(ax, moor_df, ctd_df, tide_df, mm=12, do_leg=False, C=C, M=M):
    # add tide station locations
    tide_df.plot(x='Longitude', y='Latitude',marker=M[0],ls='',mec='k',
        c=C[0], ms=mm,ax=ax,legend=False)
    # add CTD station locations
    ctd_df.plot(x='Longitude', y='Latitude',marker=M[1],ls='',mec='k',
        c=C[1], ms=mm*1.2,ax=ax,legend=False)
    # add mooring station locations
    moor_df.plot(x='Longitude', y='Latitude',marker=M[2],ls='',mec='k',
        c=C[2], ms=mm,ax=ax,legend=False)
        
    if do_leg:
        ax.plot(.95,.1,marker=M[0], c=C[0], ms=mm,mec='k', transform=ax.transAxes)
        ax.plot(.95,.075,marker=M[1], c=C[1], ms=mm*1.2,mec='k', transform=ax.transAxes)
        ax.plot(.95,.05,marker=M[2], c=C[2], ms=mm,mec='k', transform=ax.transAxes)
        ax.text(.93, .1, 'Tide Station', size=.8*fs2, ha='right',va='center',transform=ax.transAxes)
        ax.text(.93, .075, 'CTD Station', size=.8*fs2, ha='right',va='center',transform=ax.transAxes)
        ax.text(.93, .05, 'Mooring', size=.8*fs2, ha='right',va='center',transform=ax.transAxes)

# PLOTTING
plt.rc('font', size=fs)
plt.close('all')
fig = plt.figure(figsize=(16,11))

# -------------- FULL MAP -----------------------------------------------

ax = fig.add_subplot(121)
cs = ax.pcolormesh(lon[6:-6,6:],lat[6:-6,6:],v[6:-6,6:], cmap=cmap, vmin=vmin, vmax=vmax)
nudge_alpha = .1
ax.pcolormesh(lon[:,:6],lat[:,:6],v[:,:6], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[:6,:],lat[:6,:],v[:6,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
ax.pcolormesh(lon[-6:,:],lat[-6:,:],v[-6:,:], cmap=cmap, vmin=vmin, vmax=vmax, alpha=nudge_alpha)
pfun.add_bathy_contours(ax, ds, txt=False)
pfun.add_coast(ax, color='gray')
ax.axis(pfun.get_aa(ds))
pfun.dar(ax)

add_sta_loc(ax, moor_df, ctd_df, tide_df, do_leg=True)

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
# Inset colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=3)
cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
# ax.text(.08, .53, r'SST $[^\circ C]$', transform=ax.transAxes)
# ax.text(.08, .03, dstr, size=fs2, transform=ax.transAxes, style='italic')

ax.text(-123.072,46.7866,'Washington', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)
ax.text(-122.996,44.5788,'Oregon', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)
    
ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
    style='italic',ha='center',va='center',rotation=-45)
ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
    style='italic',ha='center',va='center',rotation=-10)


# -------------- FOCUS MAP -----------------------------------------------

ax = fig.add_subplot(122)
cs = ax.pcolormesh(lon,lat,v, cmap=cmap, vmin=vmin2, vmax=vmax2)
pfun.add_coast(ax, color='gray')
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

add_sta_loc(ax, moor_df, ctd_df, tide_df)

# larger station marks
mms = 20
nud = -.01
sta = 'Seattle'
x = tide_df.loc[sta,'Longitude']
y = tide_df.loc[sta,'Latitude']
ax.plot(x,y,marker=M[0],c=C[0], ms=mms,mec='k')
ax.text(x,y+nud,'T',size=fs3,ha='center',va='center', weight='bold',color='w')

sta = 'SAR003'
x = ctd_df.loc[sta,'Longitude']
y = ctd_df.loc[sta,'Latitude']
ax.plot(x,y,marker=M[1],c=C[1], ms=mms*1.5,mec='k')
ax.text(x,y+nud,'C',size=fs3,ha='center',va='center', weight='bold',color='w')

sta = 'KL027'
x = moor_df.loc[sta,'Longitude']
y = moor_df.loc[sta,'Latitude']
ax.plot(x,y,marker=M[2],c=C[2], ms=mms,mec='k')
ax.text(x,y+nud,'M',size=fs3,ha='center',va='center', weight='bold',color='w')

sta = 'PUG1718'
x = moor_df.loc[sta,'Longitude']
y = moor_df.loc[sta,'Latitude']
ax.plot(x,y,marker=M[2],c=C[2], ms=mms,mec='k')
ax.text(x,y+nud,'V',size=fs3,ha='center',va='center', weight='bold',color='w')


ax.text(.93,.97,'(b) Salish Sea', ha='right', va='top', size=fs,
    transform=ax.transAxes, weight='bold')
# Inset colorbar
cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=3)
cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
# ax.text(.08, .53, r'SST $[^\circ C]$', transform=ax.transAxes)
# ax.text(.08, .03, dstr, size=fs2, transform=ax.transAxes, style='italic')

# add labels
ax.text(-122.682,49.335,'Fraser\nRiver',size=fs2,
    style='italic',ha='center',va='center',rotation=0)
ax.text(-123.785,49.2528,'Strait of Georgia',size=fs2,
    style='italic',ha='center',va='center',rotation=-30,c='w')
ax.text(-123.434,48.2381,'Juan de Fuca',size=fs2,
    style='italic',ha='center',va='center',rotation=-0)
ax.text(-123.3,47.6143,'Puget\nSound',size=fs2,
    style='italic',ha='center',va='center',rotation=+55)
ax.text(-122.3,48.48,'Skagit\nRiver',size=fs3,
    style='italic',ha='center',va='center',
    bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
ax.text(-123.173,48.4457,'Haro\nStrait',size=fs3,
    style='italic',ha='center',va='center',
    bbox=dict(facecolor='w', edgecolor='None',alpha=.3))
    
fig.tight_layout()
fig.savefig(out_dir + 'Fig_1_stations.png')
        
ds.close()

plt.show()
plt.rcdefaults()




