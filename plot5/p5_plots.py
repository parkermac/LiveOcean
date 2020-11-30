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
    m_color = 'darkorange' # for the mooring
    
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
    pfun.add_bathy_contours(ax, ds, txt=False)
    pfun.add_coast(ax, color='gray')
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.text(.95, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
        
    ax.plot(in_dict['m_lon'], in_dict['m_lat'],'o', ms=10, c=m_color)
    
    ax.set_xticks(range(-129,-121,2))
    ax.set_yticks(range(44,52,2))
    
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['tm'])
    iot = zfun.find_nearest_ind(in_dict['ot_vec'], T['ocean_time'])
    
    ax.text(-123.08,46.64,'10 MPH\nWind', c='gray', ha='center', va='center', style='italic',
        bbox=dict(facecolor='w', edgecolor='None', alpha=.5))
    pfun.add_wind(ax, in_dict['m_lon'], in_dict['m_lat'],
        in_dict['uwind_vec'][iot], in_dict['vwind_vec'][iot])
        
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ot_vec = in_dict['ot_vec']
    zeta_vec = in_dict['zeta_vec']
    zeta_vec = zeta_vec*3.28084 # convert to feet
    
    
    ax.plot(ot_vec, zeta_vec, '-', lw=2, c=m_color)
    ax.plot(T['ocean_time'], zeta_vec[iot],'o', ms=20, c=m_color)
    x0 = ot_vec[0]; x1 = ot_vec[-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta_vec.min())-1; y1 = np.ceil(zeta_vec.max())+1
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], m_color, alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, 3,
        '+3', va='center', ha='center', c=m_color, weight='bold')
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, -3,
        '-3', va='center', ha='center', c=m_color, weight='bold')
        
    
    ax.set_yticks([])    
    ax.set_xticks([])
    
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I %p') + ' [' + dt_local.tzname() + ']',
        weight='bold', transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c=m_color,
        transform=ax.transAxes)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
        
