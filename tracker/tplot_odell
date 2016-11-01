#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:25:13 2016

Plot results of tracker.

Specifically to plot results of multiple experiments for Anthony Odell.
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matplotlib.pyplot as plt
import pickle

plp = os.path.abspath('../plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'

# make sure the output directory exists
outdir = Ldir['LOo'] + 'trackplots_odell2/'
Lfun.make_dir(outdir)

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if (m[-2:] == '.p') and ('odell2' in m):
        m_list.append(m)

for inname in m_list:
    out_fn = outdir + 'odell_5day_' + inname.split('_')[-2] + '.png'
    print('Printing ' + out_fn)
    P, G, S, PLdir = pickle.load( open( indir + inname, 'rb' ) )    
    NT, NP = P['lon'].shape    
    lonp = G['lon_psi']
    latp = G['lat_psi']
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]    
    # PLOTTING    
    plt.close()
    fig = plt.figure(figsize=(12,12))    
    ax = fig.add_subplot(111)
    pfun.add_coast(ax)    
    pfun.add_bathy_contours(ax, G, txt=True)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.text(.06, .04, ' '.join(inname.split('_')),
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical') 
    ax.set_title('Start ' + inname.split('_')[-2])
    # add the tracks
    ax.plot(P['lon'], P['lat'], '-k', alpha = 0.1)
    beach_mask = P['u'][-1,:] == 0
    ax.plot(P['lon'][:,beach_mask], P['lat'][:,beach_mask], '-r', linewidth=1)
    ax.plot(P['lon'][0,beach_mask],P['lat'][0,beach_mask],'or',
            markersize=5, alpha = .4)    
    plt.savefig(out_fn)
