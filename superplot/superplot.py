"""
Superplot development code.  Eventually this can be transfrred into a
function in plotting/roms_plots.py.

"""

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cas4', tag='v2')
Ldir['gtagex'] = Ldir['gtag'] + '_lo6biom'
import zrfun
import zfun

from importlib import reload

pth = os.path.abspath('../plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun
reload(pfun)
import pinfo

import numpy as np
import pandas as pd
import pickle
import netCDF4 as nc
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

vn = 'phytoplankton'
vlims = (0, 40)
cmap = 'jet'

# get model fields
fn = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2017.07.20/ocean_his_0001.nc'
in_dict = dict()
in_dict['fn'] = fn
ds = nc.Dataset(fn)

# get forcing fields
ffn = Ldir['LOo'] + 'superplot/forcing_' + Ldir['gtagex'] + '_2017.p'
fdf = pd.read_pickle(ffn)
fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5


# get section
G, S, T = zrfun.get_basic_info(in_dict['fn'])
# read in a section (or list of sections)
tracks_path = Ldir['data'] + 'tracks_new/'
#tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
tracks = ['Line_ps_main_v0.p']
zdeep = -300
xx = np.array([])
yy = np.array([])
for track in tracks:
    track_fn = tracks_path + track
    # get the track to interpolate onto
    pdict = pickle.load(open(track_fn, 'rb'))
    xx = np.concatenate((xx,pdict['lon_poly']))
    yy = np.concatenate((yy,pdict['lat_poly']))
for ii in range(len(xx)-1):
    x0 = xx[ii]
    x1 = xx[ii+1]
    y0 = yy[ii]
    y1 = yy[ii+1]
    nn = 20
    if ii == 0:
        x = np.linspace(x0, x1, nn)
        y = np.linspace(y0,y1, nn)
    else:
        x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
        y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

# PLOTTING

plt.close('all')
fig = plt.figure(figsize=(17,9))
fs = 18 # fontsize

# map field
ax = fig.add_subplot(131)
lon = ds['lon_psi'][:]
lat = ds['lat_psi'][:]
v =ds[vn][0, -1, 1:-1, 1:-1]
fac=pinfo.fac_dict[vn]
vv = fac * v
vv[:, :6] = np.nan
vv[:6, :] = np.nan
cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
pfun.add_coast(ax)
ax.axis(pfun.get_aa(ds))
pfun.dar(ax)
ax.set_axis_off()
# add a box for the subplot
aa = [-123.5, -122.1, 47.03, 48.8]
pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)

ax.text(.95, .07, 'LiveOcean\nPHYTOPLANKTON\n'
    + datetime.strftime(T['tm'], '%m/%d/%Y'), fontsize=fs, color='k',
    transform=ax.transAxes, horizontalalignment='center',
    fontweight='bold')

# PS map
# map field
#ax =  plt.subplot2grid((3,3), (1,2), rowspan=2)
ax = fig.add_subplot(132)
# lon = ds['lon_psi'][:]
# lat = ds['lat_psi'][:]
# v =ds[vn][0, -1, 1:-1, 1:-1]
lon = ds['lon_rho'][:]
lat = ds['lat_rho'][:]
v =ds[vn][0, -1, :, :]
fac=pinfo.fac_dict[vn]
vv = fac * v
vv[:, :6] = np.nan
vv[:6, :] = np.nan
cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1],
    cmap=cmap, shading='gouraud')
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
ax.set_axis_off()
# add section track
ax.plot(x, y, linestyle=':', color='violet', linewidth=3)


# section
ax =  fig.add_subplot(433)
#ax.plot(dist, v2['zbot'], ':k', linewidth=2)
ax.plot(dist, v2['zeta']+5, linestyle=':', color='violet', linewidth=3)
ax.set_xlim(dist.min(), dist.max())
ax.set_ylim(zdeep, 10)
sf = pinfo.fac_dict[vn] * v3['sectvarf']
# plot section
cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                   vmin=vlims[0], vmax=vlims[1], cmap=cmap)
ax.set_axis_off()
ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
    transform=ax.transAxes)

fig.tight_layout()
plt.show()

# get the day
tm = T['tm'] # datetime
TM = datetime(tm.year, tm.month, tm.day)

# get yearday
yearday = fdf['yearday'].values
this_yd = fdf.loc[TM, 'yearday']

# Tides
ax = fig.add_subplot(436)
ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
    lw=3, alpha=.4)
ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
    marker='o', color='r', markersize=7)
ax.set_xlim(0,365)
ax.set_ylim(.4, 1.7)
ax.set_axis_off()
ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
    alpha=.4, fontsize=fs, horizontalalignment='right')
ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
    alpha=.4, fontsize=fs, horizontalalignment='right')


# wind
alpha=.5
ax = fig.add_subplot(439)
w = fdf['8-day NS Wind Stress (Pa)'].values
wp = w.copy()
wp[w<0] = np.nan
wm = w.copy()
wm[w>0] = np.nan
tt = np.arange(len(w))
ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
    marker='o', color='r', markersize=7)
ax.set_axis_off()
ax.set_xlim(0,365)
ax.set_ylim(-.15, .25)
ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
    color='g', alpha=alpha, fontsize=fs)
ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
    color='b', alpha=alpha, fontsize=fs)


# Rivers
cr = fdf['Columbia R. Flow (1000 m3/s)'].values
fr = fdf['Fraser R. Flow (1000 m3/s)'].values
sr = fdf['Skagit R. Flow (1000 m3/s)'].values
this_yd = fdf.loc[TM, 'yearday']
ax = fig.add_subplot(4,3,12)
# lw = 3
# ax.plot(yearday, cr, color='orange', lw=lw)
# ax.plot(yearday, fr, color='violet', lw=lw)
# ax.plot(yearday, sr, color='brown', lw=lw)
alpha = .6
ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)

ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
    marker='o', color='r', markersize=7)
ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
    marker='o', color='r', markersize=7)
ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
    marker='o', color='r', markersize=7)
    
ax.text(.8, .85, 'Columbia River', transform=ax.transAxes,
    color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
ax.text(.8, .70, 'Fraser River', transform=ax.transAxes,
    color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
ax.text(.8, .55, 'Skagit River', transform=ax.transAxes,
    color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)

ax.set_axis_off()
    
#ax.set_xlabel('Yearday 2017')
ax.set_xlim(0,365)    
ax.set_ylim(0,20)

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')


