"""
Module of plotting functions.
"""

# imports

# The calling function, p5.py, has already put alpha on the path.
import Lfun
Ldir = Lfun.Lstart()
if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import zfun
import zrfun

from importlib import reload
import pfun
reload(pfun)
import pinfo
reload(pinfo)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def P_1(Q, M):
    fig = plt.figure(figsize=(6.5,12))
    fs=18
    plt.rc('font', size=fs)
    
    ds = nc.Dataset(Q['fn'])
    if Q['dom'] == 'full':
        Q['aa'] = pfun.get_aa(ds)
    aa = Q['aa']
    
    vn = Q['vn']
    if Q['bot']:
        nlev = 0
        tstr = 'Bottom'
    else:
        nlev = -1
        tstr = 'Surface'
        
    if vn == 'speed':
        u = ds['u'][0,nlev,:,:]
        v = ds['v'][0,nlev,:,:]
        u[u.mask] = 0
        v[v.mask] = 0
        fld0 = 0 * ds['salt'][0,-1,1:-1,1:-1]
        uu = (u[1:-1,:-1] + u[1:-1,1:])/2
        vv = (v[:-1,1:-1] + v[1:,1:-1])/2
        fld = np.sqrt(uu*uu + vv*vv)
        fld = np.ma.masked_where(fld0.mask, fld)
    else:
        fld = ds[vn][0,nlev,1:-1,1:-1]*pinfo.fac_dict[vn]
        
    if Q['emask']:
        fld = pfun.mask_edges(ds, fld, Q)
    if Q['avl']:
        # set vmax and vmin if needed
        pfun.get_vlims(ds, fld, Q)
        
    # MAP FIELD
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], fld,
        cmap=pinfo.cmap_dict[vn], vmin=Q['vmin'], vmax=Q['vmax'])
    ax.text(.95, .99,'%s %s %s' % (tstr, pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    
    if Q['dom'] == 'full':
        pfun.add_bathy_contours(ax, ds, txt=False, depth_levs=[200])
        
    if vn == 'speed':
        pfun.add_velocity_vectors(ax, aa, ds, Q['fn'], v_scl=Q['v_scl'])
        
    pfun.add_coast(ax, color='k')
    ax.axis(aa)
    pfun.dar(ax)
        
    # Inset colorbar
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    
    # Mooring location
    ax.plot(M['lon'], M['lat'],'ok', ms=5)
    
    # Wind vector
    T = zrfun.get_basic_info(Q['fn'], only_T=True)
    pfun.add_wind(ax, M, T)
    pfun.add_wind_text(ax, aa, M, fs)
        
    # axes labeling
    ax.set_xticks(Q['xtl'])
    ax.set_yticks(Q['ytl'])
    ax.tick_params(labelsize=.7*fs)
    # ax.tick_params(axis="y",direction="in", pad=-28)#, labelcolor='gray')
    # ax.tick_params(axis="x",direction="in", pad=-21)#, labelcolor='gray')
    
        
    # MOORING TIME SERIES
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    pfun.plot_time_series(ax, M, T)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(Q['fn_out']) > 0:
        plt.savefig(Q['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
    
