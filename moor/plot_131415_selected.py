# -*- coding: utf-8 -*-

"""
Plots mooring records.  Created for the Pacific Anomalies Workshop 2016_01.

Inclues analysis of stratification.
"""

# setup
from datetime import datetime
import matplotlib.dates as mdates
import numpy as np

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)

from importlib import reload
import Lfun
reload(Lfun)

import zfun  # plotting functions
reload(zfun)

import seawater as sw

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'

# choose the data to plot
if False:
    print('\n%s\n' % '** Choose mooring file to plot **')
    m_list_raw = os.listdir(indir)
    m_list = []
    for m in m_list_raw:
        if m[-2:] == '.p':
            m_list.append(m)
    Npt = len(m_list)
    m_dict = dict(zip(range(Npt), m_list))
    for npt in range(Npt):
        print(str(npt) + ': ' + m_list[npt])
    my_npt = int(input('-- Input number -- '))
    inname = m_dict[my_npt]
else:
    inname = 'cascadia1_base_lo1_47N124.5W_low_pass_2013.01.02_2015.12.31.p'
  
import pickle
V, v1_list, v2_list, v3_list, G, S, Lon, Lat, sta_name, h = pickle.load( open( indir + inname, 'rb' ) )

#%% messing with z
z_rho, z_w = zfun.get_z(h, 0*h, S)

z_rho_ext = zfun.make_full((-h, z_rho, np.array([0.,])))

zees = np.array([-70, -35, -20, -10, -0])
zcolor = ['k', 'b', 'g', 'orange', 'r']
nz = len(zees)

int_arr = zfun.get_interpolant_fast(zees, z_rho_ext)

v3_list_alt = ['temp', 'salt', 'v']
v3_names = ['Temperature $(^{\circ}C)$', 'Salinity',
            'NS Velocity $(m s^{-1})$']
v3_lims = [(6,19), (34,26), (-.6,.6)]
v3_limdict = dict(zip(v3_list_alt, v3_lims))
v3_namedict = dict(zip(v3_list_alt, v3_names))

v2_list_alt = ['pe_mix_v', 'svstr', 'shflux']
v2_names = ['Energy to Mix $(Jm^{-3})$', 'NS wind Stress (Pa)',
            'Surface Net Heat Flux $(Wm^{-2})$']
v2_lims = [(0,160), (-.15,.25), (-300,300)]
v2_limdict = dict(zip(v2_list_alt, v2_lims))
v2_namedict = dict(zip(v2_list_alt, v2_names))

# this creates arrays with time series interpolated to the selected depths
VV = dict()
for vv in v3_list_alt:
    fld_ext = zfun.make_full((V[vv],))    
    VV[vv] = (fld_ext[int_arr[:,0].astype(int),:]*(1-int_arr[:,2]).reshape(nz,1) +
                fld_ext[int_arr[:,1].astype(int),:]*int_arr[:,2].reshape(nz,1))
                
#%% stratification           
sigma = sw.dens(V['salt'], V['temp'], 0.) - 1000.
N, NT = sigma.shape
dz = np.diff(z_w)
dz1 = dz.reshape(N,1)
H = dz.sum()
sig_bar = (sigma*dz1).sum(axis=0)/H
sig_prime = sigma - sig_bar.reshape(1,NT)*np.ones((N,1))
g = 9.8

# energy (J m-3) per unit volume to vertically mix the water column
V['pe_mix_v'] = -(g*sig_prime*z_rho.reshape(N,1)).sum(axis=0)/H

                
#%% time fields
                
mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
year = 2013 + (mdays - mdates.date2num(datetime(2013,1,1)))/356

#%% plotting
import matplotlib.pyplot as plt
plt.close()
NR = 2; NC = 3
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,10), squeeze=False)

cc = 0
nfilt = 20
for vn in v3_list_alt:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    for iz in range(nz):
        ax.plot(year, zfun.filt_hanning(VV[vn][iz,:],n=nfilt), '-',
                color=zcolor[iz])
    ax.set_xlim(2013, 2016)
    ax.ticklabel_format(useOffset=False, axis='x')
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_xticks([2013.5, 2014.5, 2015.5])
    ax.set_xticklabels([2013, 2014, 2015])
    ax.set_title(v3_namedict[vn])
    if vn == 'temp':
        tx0 = inname.find('47N')
        tx1 = inname.find('W')
        tx_str = inname[tx0:tx1+1]
        ax.text(.05, .9,
                tx_str,
                transform=ax.transAxes,
                fontsize=14)
        # and add the depth averaged temperature
        t_bar = (V[vn]*dz1).sum(axis=0)/H
        ax.plot(year, zfun.filt_hanning(t_bar,n=nfilt),
                '-m', linewidth=3, alpha=.5)
        ax.text(.95, .9,
                'DEPTH MEAN',
                transform=ax.transAxes,
                horizontalalignment='right',
                fontsize=14,
                fontweight='bold',
                color='m',
                alpha=.5)
    if vn == 'salt':
        ax.invert_yaxis()
        zcount = 0
        for zz in zees:
            ax.text(.95, .5 + .09*zcount,
                    str(zz) + ' m',
                    color=zcolor[zcount],
                    horizontalalignment='right',
                    transform=ax.transAxes,
                    fontsize=12)
            zcount += 1
    ax.set_ylim(v3_limdict[vn])
    ax.grid()
    yl = ax.get_ylim()
    ax.plot([2014, 2014],yl,'-k')
    ax.plot([2015, 2015],yl,'-k')
    ax.set_ylim(yl)
    if vn in ['v']:
        xl = ax.get_xlim()
        ax.plot(xl,[0,0],'-k')
    cc += 1
for vn in v2_list_alt:
    ir = int(np.floor(cc/NC))
    ic = int(cc - NC*ir)
    ax = axes[ir, ic]
    ax.plot(year, zfun.filt_hanning(V[vn],n=nfilt))
    ax.set_xlim(2013, 2016)
    ax.ticklabel_format(useOffset=False, axis='x')
    ax.ticklabel_format(useOffset=False, axis='y')
    ax.set_xticks([2013.5, 2014.5, 2015.5])
    ax.set_xticklabels([2013, 2014, 2015])
    ax.set_title(v2_namedict[vn])
    ax.set_ylim(v2_limdict[vn])
    ax.grid()
    if vn in ['svstr', 'shflux']:
        xl = ax.get_xlim()
        ax.plot(xl,[0,0],'-k')
    yl = ax.get_ylim()
    ax.plot([2014, 2014],yl,'-k')
    ax.plot([2015, 2015],yl,'-k')
    ax.set_ylim(yl)
    cc += 1
plt.show()

#%% printing
plt.savefig(Ldir['LOo']+ 'moor/' + tx_str + '.png')
