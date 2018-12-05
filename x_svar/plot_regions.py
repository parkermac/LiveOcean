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

pdir = Ldir['parent'] + 'ptools_output/polygons/'

p_list = ['basin_hood_canal','basin_main_basin',
    'basin_south_sound','basin_whidbey',
    'strait_admiralty','strait_tacoma_narrows',
    'strait_hood_canal','strait_triple_junction']

v_dict = dict()
for pn in p_list:
    poly_fn = pdir + pn + '.p'
    poly = pickle.load(open(poly_fn, 'rb'))
    px = poly['lon_poly']
    py = poly['lat_poly']
    v = np.ones((len(px),2))
    v[:,0] = px
    v[:,1] = py
    v_dict[pn] = v

# PLOTTING
plt.close('all')


# MAPS OF NET MIXING IN REGIONS
fig = plt.figure(figsize=(8,10))
aa = [-123.25, -122, 47, 48.5]

ax = fig.add_subplot(111)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Puget Sound Basins & Straits')

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

fig.tight_layout()


plt.show()


