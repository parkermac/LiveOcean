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
import p5_fun as pfun
reload(pfun)
import p5_info as pinfo
reload(pinfo)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pickle
from datetime import datetime, timedelta
import pandas as pd

def P_salt_1(in_dict):
    # single panel of salinity

    # START
    fig = plt.figure(figsize=(6.5,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    fs=18
    plt.rc('font', size=fs)
    vn = 'salt'
    
    if True:
        import cmocean
        cmap = cmocean.cm.haline
    else:
        cmap = 'nipy_spectral'
    
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    vlims_fac=.5
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=cmap, fac=pinfo.fac_dict[vn], vlims_fac=vlims_fac)
    #cb = fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=False)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.text(.9, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    
    #pfun.add_info(ax, in_dict['fn'], fs=fs)
    
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    
    ax.set_xticks(range(-129,-121,2))
    ax.set_yticks(range(44,52,2))
    
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ax.plot(in_dict['ot_vec'], in_dict['zeta_vec'])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
        
