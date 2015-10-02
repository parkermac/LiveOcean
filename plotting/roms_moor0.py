"""
Extract and plot a mooring-like record..
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta

import zfun; reload(zfun) # plotting functions
import matfun; reload(matfun) # functions for working with mat files

which_home = os.environ.get("HOME") # This works even when called by cron.
if which_home == '/Users/PM5': # mac version
    do_mac = True
elif which_home == '/home/parker': # fjord version
    do_mac = False
    
def make_date_list(dt0,dt1,Ldir): # a helpful function
    del_dt = timedelta(1)
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + del_dt
    return date_list

list_type = 'hindcast'
fn_list = []
if list_type == 'hindcast':
    dt0 = datetime(2013,6,1) # first day
    dt1 = datetime(2013,6,30) # last day
    Ldir['gtag'] = Ldir['gtag'] + '_2013'
    date_list = make_date_list(dt0,dt1,Ldir)
    
# target position (-124.5, 47 = RN)
Lon = np.array([-124.5])
Lat = np.array([47.])

indir = Ldir['roms'] + 'output/' + Ldir['gtag'] + '/f' + date_list[0] + '/'
fn = indir + 'ocean_his_0002.nc'
[G] = zfun.get_basic_info(fn, getS=False, getT=False)

# get interpolants for this point
Xit = dict(); Yit = dict()
Aix = dict(); Aiy = dict()
for grd in ['rho', 'u', 'v']:
    xx = G['lon_' + grd][1,:]
    yy = G['lat_' + grd][:,1]
    xit_list = zfun.get_interpolant(Lon, xx)
    yit_list = zfun.get_interpolant(Lat, yy)
    # this just pulls the interpolant tuple out of the (one-element) list
    # that get_interpolant returns
    Xit[grd] = xit_list[0]; Yit[grd] = yit_list[0]
    # create little arrays that are used in the actual interpolation
    Aix[grd] = np.array([1-Xit[grd][2], Xit[grd][2]]).reshape((1,1,2))           
    Aiy[grd] = np.array([1-Yit[grd][2], Yit[grd][2]]).reshape((1,2))

v1_list = ['ocean_time']
v2_list = ['zeta','sustr','svstr','swrad','lwrad',
    'shflux','latent','sensible']
v3_list = ['temp','salt','u','v']

V = dict()
for vv in v1_list:
    V[vv] = np.array([])
for vv in v2_list:
    V[vv] = np.array([])
for vv in v3_list:
    V[vv] = np.array([])
    
def get_its(ds, vv, Xit, Yit, Aix, Aiy):
    dims = ds.variables[vv].dimensions
    if 'eta_rho' in dims:
        grd = 'rho'
    elif 'eta_u' in dims:
        grd = 'u'
    elif 'eta_v' in dims:
        grd = 'v'
    else:
        print 'grid error!'
    xit = Xit[grd]; yit = Yit[grd]
    aix = Aix[grd]; aiy = Aiy[grd]
    return xit, yit, aix, aiy

count = 0    
for dd in date_list:
    print 'Working on date_list item: ' + dd
    sys.stdout.flush()
    if list_type == 'hindcast':
        indir = Ldir['roms'] + 'output/' + Ldir['gtag'] + '/f' + dd + '/'
        ds = nc.MFDataset(indir + 'ocean_his*.nc')
    elif list_type == 'pnwtox':
            ds = nc.MFDataset(indir + 'ocean_his_' + dd + '*.nc')
    for vv in v1_list:
        vtemp = ds.variables[vv][:].squeeze()
        V[vv] = np.append(V[vv], vtemp)
    for vv in v2_list:
        xit, yit, aix, aiy = get_its(ds, vv, Xit, Yit, Aix, Aiy)
        vvtemp = ds.variables[vv][:, yit[:2], xit[:2]].squeeze()
        vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
        V[vv] = np.append(V[vv], vtemp)        
    for vv in v3_list:
        xit, yit, aix, aiy = get_its(ds, vv, Xit, Yit, Aix, Aiy)
        vvtemp = ds.variables[vv][:, -1, yit[:2], xit[:2]].squeeze()
        vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
        V[vv] = np.append(V[vv], vtemp)
    # listing of contents, if needed
    if count == 0 and True:   
        zfun.ncd(ds)
    count += 1
    ds.close()
                   
# plotting
import matplotlib.pyplot as plt

plt.close()

fig = plt.figure(figsize=(14, 10))

days = (V['ocean_time'] - V['ocean_time'][0])/86400.

vn_list = ['zeta','temp','salt',
    'v','swrad','lwrad',
    'shflux','latent','sensible']

axn = 1
for vn in vn_list:    
    ax = fig.add_subplot(3,3,axn)
    ax.plot(days, V[vn])
    ax.set_xlim(0, days.max())
    ax.set_title(vn)
    axn += 1

plt.show()

#plt.savefig(Ldir['LOo'] + 'plots/' + list_type + '.png')