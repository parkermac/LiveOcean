"""
Plots layer records for svar extractions.
"""

# setup
import matplotlib.pyplot as plt
import numpy as np
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

indir = Ldir['LOo'] + 'svar/'
in_fn = indir + 'spring_neap_2017.p'
dd = pickle.load(open(in_fn, 'rb'))

xp = dd['xp']
yp = dd['yp']
x = dd['x']
y = dd['y']
aa = dd['aa']
DA = dd['DA']
mmix_neap = dd['mmix_neap']
mmix_spring = dd['mmix_spring']

mmix = mmix_neap + mmix_spring

mmixda_neap = mmix_neap * DA
mmixda_spring = mmix_spring * DA

xy2 = np.stack((x.flatten(),y.flatten()), axis=1)
M_neap = mmixda_neap.flatten()
M_spring = mmixda_spring.flatten()

# form some area integrals
import matplotlib.path as mpath
v_dict = dict() # V is for the Vertices of the polygons
p_dict = dict()
M_neap_dict = dict()
M_spring_dict = dict()
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
    M_neap_dict[pn] = M_neap[p_in].sum()
    M_spring_dict[pn] = M_spring[p_in].sum()

lmix = np.log10(mmix+1e-8)
lmix_neap = np.log10(mmix_neap+1e-8)
lmix_spring = np.log10(mmix_spring+1e-8)

# PLOTTING
plt.close('all')

# MAPS OF MIXING per unit area
fig = plt.figure(figsize=(18,8))
vmin = -6; vmax = -2

# full region
ax = fig.add_subplot(131)
cs = ax.pcolormesh(xp, yp, lmix[1:-1, 1:-1], vmin=vmin, vmax=vmax, cmap='rainbow')
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_title('$log_{10}[Average\ Mixing\ (psu^{2}\ m\ s^{-1})]$')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# close up of Puget Sound
aa = [-123.25, -122, 47, 48.5]
ax = fig.add_subplot(132)
cs = ax.pcolormesh(xp, yp, lmix_neap[1:-1, 1:-1], vmin=vmin, vmax=vmax, cmap='rainbow')
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_xlabel('Longitude')
ax.set_title('Close-up of Puget Sound: NEAP')

ax = fig.add_subplot(133)
cs = ax.pcolormesh(xp, yp, lmix_spring[1:-1, 1:-1], vmin=vmin, vmax=vmax, cmap='rainbow')
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis(aa)
ax.set_xlabel('Longitude')
ax.set_title('Close-up of Puget Sound: SPRING')



fig.tight_layout()

# MAPS OF NET MIXING IN REGIONS
fig = plt.figure(figsize=(15,8))
aa = [-123.25, -122, 47, 48.5]

ax = fig.add_subplot(131)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Basins & Straits')

# add region polygons and names
clr_dict = dict()
for pn in p_list:
    v = v_dict[pn]
    xpoly = v[:,0]; ypoly = v[:,1]
    xpoly = np.concatenate((xpoly, np.array(xpoly[0]).reshape(1,)))
    ypoly = np.concatenate((ypoly, np.array(ypoly[0]).reshape(1,)))
    lh = ax.plot(xpoly,ypoly, '-', linewidth=3, alpha=.5)
    clr_dict[pn] = lh[0].get_color()
    pstr = pn.split('_')
    pstr = ' '.join(pstr[1:]).title()
    ax.text(xpoly.mean(), ypoly.mean(), pstr,
        fontweight='bold', horizontalalignment='center',
        verticalalignment='center', color=clr_dict[pn], fontsize=14)

# close ups of Puget Sound

ax = fig.add_subplot(132)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_title('Net Mixing NEAP $(10^{3}\ psu^{2}\ m^{3}\ s^{-1})$')
# add mixing info
for pn in p_list:
    v = v_dict[pn]
    xpoly = v[:,0]; ypoly = v[:,1]
    xpoly = np.concatenate((xpoly, np.array(xpoly[0]).reshape(1,)))
    ypoly = np.concatenate((ypoly, np.array(ypoly[0]).reshape(1,)))
    mstr = '%d' % (int(M_neap_dict[pn]/1000))
    ax.text(xpoly.mean(), ypoly.mean(), mstr,
        fontweight='bold', horizontalalignment='center',
        verticalalignment='center', color=clr_dict[pn], fontsize=28)
        
ax = fig.add_subplot(133)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_title('Net Mixing SPRING $(10^{3}\ psu^{2}\ m^{3}\ s^{-1})$')
# add mixing info
for pn in p_list:
    v = v_dict[pn]
    xpoly = v[:,0]; ypoly = v[:,1]
    xpoly = np.concatenate((xpoly, np.array(xpoly[0]).reshape(1,)))
    ypoly = np.concatenate((ypoly, np.array(ypoly[0]).reshape(1,)))
    mstr = '%d' % (int(M_spring_dict[pn]/1000))
    ax.text(xpoly.mean(), ypoly.mean(), mstr,
        fontweight='bold', horizontalalignment='center',
        verticalalignment='center', color=clr_dict[pn], fontsize=28)

fig.tight_layout()


plt.show()


