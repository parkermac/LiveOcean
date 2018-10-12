"""
Plots layer records for svar extractions.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir0 = Ldir['LOo'] + 'layer/'

# choose the layer extraction to plot
if False:
    item = Lfun.choose_item(indir0)
else:
    item = 'cas4_v2_lo6biom_2017.07.20_2017.07.22'
indir = indir0 + item + '/'
if False:
    infile = Lfun.choose_item(indir, tag='.nc')
    fn = indir + infile
else:
    fn = indir + 'svar_hourly.nc'

ds = nc.Dataset(fn)
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]
aa = pfun.get_aa(ds)
mix = ds['mix'][:]
mask= ds['mask_rho'][:] # 1 on water, 0 on land
DA = ds['DA'][:] # cell horizontal area
x = ds['lon_rho'][:]
y = ds['lat_rho'][:]
ds.close()

days = (ot - ot[0])/86400.
mdays = Lfun.modtime_to_mdate_vec(ot)
mdt = mdates.num2date(mdays) # list of datetimes of data

mmix = mix.mean(axis=0)

mmixda = mmix * DA


xy2 = np.stack((x.flatten(),y.flatten()), axis=1)
M = mmixda.flatten()

# form some area integrals
import matplotlib.path as mpath
v_dict = dict() # V is for the Vertices of the polygons
p_dict = dict()
M_dict = dict()
pdir = Ldir['parent'] + 'ptools_output/polygons/'

p_list = ['basin_hood_canal','basin_main_basin',
    'basin_south_sound','basin_whidbey',
    'strait_admiralty','strait_tacoma_narrows',
    'strait_hood_canal','strait_triple_junction']

for pn in p_list:
    poly_fn = pdir + pn + '.p'
    poly = pickle.load(open(poly_fn, 'rb'))
    px = poly['lon_poly']
    py = poly['lat_poly']
    v = np.ones((len(px),2))
    v[:,0] = px
    v[:,1] = py
    v_dict[pn] = v
    p_dict[pn] = mpath.Path(v)
for pn in p_list:
    p = p_dict[pn]
    p_in = p.contains_points(xy2) # boolean
    M_dict[pn] = M[p_in].sum()

lmix = np.log10(mmix+1e-8)

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(11,8))

vmin = -6; vmax = -3

# full region
ax = fig.add_subplot(121)
cs = ax.pcolormesh(xp, yp, lmix[1:-1, 1:-1], vmin=vmin, vmax=vmax, cmap='rainbow')
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('$log_{10}(Mixing)\ (W\ m^{-2})$')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# close up of Puget Sound
aa = [-123.5, -122, 47, 49]
ax = fig.add_subplot(122)
cs = ax.pcolormesh(xp, yp, lmix[1:-1, 1:-1], vmin=vmin, vmax=vmax, cmap='rainbow')
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_xlabel('Longitude')
ax.set_title('Close-up of Puget Sound')

# add mixing info
for pn in p_list:
    v = v_dict[pn]
    xp = v[:,0]; yp = v[:,1]
    xp = np.concatenate((xp, np.array(xp[0]).reshape(1,)))
    yp = np.concatenate((yp, np.array(yp[0]).reshape(1,)))
    lh = ax.plot(xp,yp, '-', linewidth=3)
    clr = lh[0].get_color()
    ax.text(xp.mean(), yp.mean(), pn,
        fontweight='bold', horizontalalignment='center', color=clr)
    ax.text(xp.mean(), yp.mean()-.05,
        '%0.1f kilowatts' % (M_dict[pn]/1000),
            fontweight='bold', horizontalalignment='center', color='k')

fig.tight_layout()

plt.show()


