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
from datetime import datetime, timedelta
import pandas as pd
import cmocean

def P_full_salt(in_dict):
    fig = plt.figure(figsize=(6.5,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    fs=18
    plt.rc('font', size=fs)
    vn = 'salt'
    vlims_fac=.5
    cmap = cmocean.cm.haline
    
    if len(in_dict['aa']) == 0:
        aa = pfun.get_aa(ds)
    else:
        aa = in_dict['aa']
    
    fld = ds[vn][0,-1,1:-1,1:-1]
    fld = fld * pinfo.fac_dict[vn]
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    if in_dict['auto_vlims']:
        (vmin, vmax) = pfun.auto_lims(fld, vlims_fac=vlims_fac)
        in_dict['vmin'] = vmin
        in_dict['vmax'] = vmax
    else:
        vmin = in_dict['vmin']
        vmax = in_dict['vmax']
    
    # plot the map field
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = ax.pcolormesh(x, y, fld, cmap=cmap, vmin=vmin, vmax=vmax)
    pfun.add_bathy_contours(ax, ds, txt=False, depth_levs=[200])
    pfun.add_coast(ax, color='gray')
    ax.axis(aa)
    pfun.dar(ax)
    ax.text(.95, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    # mooring location
    mcolor = 'firebrick' # color for the mooring
    ax.plot(in_dict['m_lon'], in_dict['m_lat'],'o', ms=5, c=mcolor)
    # axes labeling
    ax.set_xticks(in_dict['xtr'])
    ax.set_yticks(in_dict['ytr'])
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    # wind speed
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['tm'])
    iot = zfun.find_nearest_ind(in_dict['ot_vec'], T['ocean_time'])
    wscl = 20
    wcolor = 'purple'
    ax.text(in_dict['xwt'], in_dict['ywt'],str(wscl)+' knot\nwind', c=wcolor,
        ha='center', va='center', style='italic',
        bbox=dict(facecolor='w', edgecolor='None', alpha=.8))
    pfun.add_wind(ax, in_dict['m_lon'], in_dict['m_lat'],
        in_dict['uwind_vec'][iot], in_dict['vwind_vec'][iot], scl=wscl, color=wcolor)
        
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ot_vec = in_dict['ot_vec']
    zeta_vec = in_dict['zeta_vec']
    zeta_vec = zeta_vec*3.28084 # convert meters to feet
    ax.plot(ot_vec, zeta_vec, '-', lw=2, c=mcolor)
    ax.plot(T['ocean_time'], zeta_vec[iot],'o', ms=20, c=mcolor)
    x0 = ot_vec[0]; x1 = ot_vec[-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta_vec.min())-5
    y1 = np.ceil(zeta_vec.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], mcolor, alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p') + ' ' + dt_local.tzname() + '',
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c=mcolor,
        transform=ax.transAxes)
    # add day/night
    Srise = in_dict['Srise']
    Sset = in_dict['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
    
def P_PS_salt(in_dict):
    fig = plt.figure(figsize=(6.5,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    fs=18
    plt.rc('font', size=fs)
    vn = 'salt'
    vlims_fac=.5
    cmap = cmocean.cm.haline
    
    if len(in_dict['aa']) == 0:
        aa = pfun.get_aa(ds)
    else:
        aa = in_dict['aa']
    
    fld = ds[vn][0,-1,1:-1,1:-1]
    fld = fld * pinfo.fac_dict[vn]
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    if in_dict['auto_vlims']:
        (vmin, vmax) = pfun.auto_lims(fld, vlims_fac=vlims_fac)
        in_dict['vmin'] = vmin
        in_dict['vmax'] = vmax
    else:
        vmin = in_dict['vmin']
        vmax = in_dict['vmax']
    
    # plot the map field
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = ax.pcolormesh(x, y, fld, cmap=cmap, vmin=vmin, vmax=vmax)
    #pfun.add_bathy_contours(ax, ds, txt=False, depth_levs=[200])
    pfun.add_coast(ax, color='gray')
    ax.axis(aa)
    pfun.dar(ax)
    ax.text(.95, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    # mooring location
    mcolor = 'firebrick' # color for the mooring
    ax.plot(in_dict['m_lon'], in_dict['m_lat'],'o', ms=5, c=mcolor)
    # axes labeling
    ax.set_xticks(in_dict['xtr'])
    ax.set_yticks(in_dict['ytr'])
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    # wind speed
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['tm'])
    iot = zfun.find_nearest_ind(in_dict['ot_vec'], T['ocean_time'])
    wscl = 20
    wcolor = 'purple'
    ax.text(in_dict['xwt'], in_dict['ywt'],str(wscl)+' knot\nwind', c=wcolor,
        ha='center', va='center', style='italic',
        bbox=dict(facecolor='w', edgecolor='None', alpha=.8))
    pfun.add_wind(ax, in_dict['m_lon'], in_dict['m_lat'],
        in_dict['uwind_vec'][iot], in_dict['vwind_vec'][iot], scl=wscl, color=wcolor)
        
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ot_vec = in_dict['ot_vec']
    zeta_vec = in_dict['zeta_vec']
    zeta_vec = zeta_vec*3.28084 # convert meters to feet
    ax.plot(ot_vec, zeta_vec, '-', lw=2, c=mcolor)
    ax.plot(T['ocean_time'], zeta_vec[iot],'o', ms=20, c=mcolor)
    x0 = ot_vec[0]; x1 = ot_vec[-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta_vec.min())-5
    y1 = np.ceil(zeta_vec.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], mcolor, alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p') + ' ' + dt_local.tzname() + '',
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c=mcolor,
        transform=ax.transAxes)
    # add day/night
    Srise = in_dict['Srise']
    Sset = in_dict['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
    
def P_full_temp(in_dict):
    fig = plt.figure(figsize=(6.5,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    fs=18
    plt.rc('font', size=fs)
    vn = 'temp'
    vlims_fac=3
    cmap = cmocean.cm.thermal
    
    if len(in_dict['aa']) == 0:
        aa = pfun.get_aa(ds)
    else:
        aa = in_dict['aa']
    
    fld = ds[vn][0,-1,1:-1,1:-1]
    fld = fld * pinfo.fac_dict[vn]
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    if in_dict['auto_vlims']:
        (vmin, vmax) = pfun.auto_lims(fld, vlims_fac=vlims_fac)
        in_dict['vmin'] = vmin
        in_dict['vmax'] = vmax
    else:
        vmin = in_dict['vmin']
        vmax = in_dict['vmax']
    
    # plot the map field
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = ax.pcolormesh(x, y, fld, cmap=cmap, vmin=vmin, vmax=vmax)
    pfun.add_bathy_contours(ax, ds, txt=False, depth_levs=[200])
    pfun.add_coast(ax, color='gray')
    ax.axis(aa)
    pfun.dar(ax)
    ax.text(.95, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    # mooring location
    mcolor = 'firebrick' # color for the mooring
    ax.plot(in_dict['m_lon'], in_dict['m_lat'],'o', ms=5, c=mcolor)
    # axes labeling
    ax.set_xticks(in_dict['xtr'])
    ax.set_yticks(in_dict['ytr'])
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    # wind speed
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['tm'])
    iot = zfun.find_nearest_ind(in_dict['ot_vec'], T['ocean_time'])
    wscl = 20
    wcolor = 'purple'
    ax.text(in_dict['xwt'], in_dict['ywt'],str(wscl)+' knot\nwind', c=wcolor,
        ha='center', va='center', style='italic',
        bbox=dict(facecolor='w', edgecolor='None', alpha=.8))
    pfun.add_wind(ax, in_dict['m_lon'], in_dict['m_lat'],
        in_dict['uwind_vec'][iot], in_dict['vwind_vec'][iot], scl=wscl, color=wcolor)
        
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ot_vec = in_dict['ot_vec']
    zeta_vec = in_dict['zeta_vec']
    zeta_vec = zeta_vec*3.28084 # convert meters to feet
    ax.plot(ot_vec, zeta_vec, '-', lw=2, c=mcolor)
    ax.plot(T['ocean_time'], zeta_vec[iot],'o', ms=20, c=mcolor)
    x0 = ot_vec[0]; x1 = ot_vec[-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta_vec.min())-5
    y1 = np.ceil(zeta_vec.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], mcolor, alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p') + ' ' + dt_local.tzname() + '',
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c=mcolor,
        transform=ax.transAxes)
    # add day/night
    Srise = in_dict['Srise']
    Sset = in_dict['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
        
def P_PS_temp(in_dict):
    fig = plt.figure(figsize=(6.5,12))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    fs=18
    plt.rc('font', size=fs)
    vn = 'temp'
    vlims_fac=3
    cmap = cmocean.cm.thermal
    
    if len(in_dict['aa']) == 0:
        aa = pfun.get_aa(ds)
    else:
        aa = in_dict['aa']
    
    fld = ds[vn][0,-1,1:-1,1:-1]
    fld = fld * pinfo.fac_dict[vn]
    x = ds['lon_psi'][:]
    y = ds['lat_psi'][:]
    
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    if in_dict['auto_vlims']:
        (vmin, vmax) = pfun.auto_lims(fld, vlims_fac=vlims_fac)
        in_dict['vmin'] = vmin
        in_dict['vmax'] = vmax
    else:
        vmin = in_dict['vmin']
        vmax = in_dict['vmax']
    
    # plot the map field
    ax = plt.subplot2grid((7,1), (0,0), rowspan=6)
    cs = ax.pcolormesh(x, y, fld, cmap=cmap, vmin=vmin, vmax=vmax)
    #pfun.add_bathy_contours(ax, ds, txt=False, depth_levs=[200])
    pfun.add_coast(ax, color='gray')
    ax.axis(aa)
    pfun.dar(ax)
    ax.text(.95, .99,'Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        transform=ax.transAxes, ha='right', va='top', weight='bold')
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    # mooring location
    mcolor = 'firebrick' # color for the mooring
    ax.plot(in_dict['m_lon'], in_dict['m_lat'],'o', ms=5, c=mcolor)
    # axes labeling
    ax.set_xticks(in_dict['xtr'])
    ax.set_yticks(in_dict['ytr'])
    ax.tick_params(axis="y",direction="in", pad=-28, labelcolor='gray')
    ax.tick_params(axis="x",direction="in", pad=-21, labelcolor='gray')
    # wind speed
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['tm'])
    iot = zfun.find_nearest_ind(in_dict['ot_vec'], T['ocean_time'])
    wscl = 20
    wcolor = 'purple'
    ax.text(in_dict['xwt'], in_dict['ywt'],str(wscl)+' knot\nwind', c=wcolor,
        ha='center', va='center', style='italic',
        bbox=dict(facecolor='w', edgecolor='None', alpha=.8))
    pfun.add_wind(ax, in_dict['m_lon'], in_dict['m_lat'],
        in_dict['uwind_vec'][iot], in_dict['vwind_vec'][iot], scl=wscl, color=wcolor)
        
    # mooring time series
    ax = plt.subplot2grid((7,1), (6,0), rowspan=1)
    ot_vec = in_dict['ot_vec']
    zeta_vec = in_dict['zeta_vec']
    zeta_vec = zeta_vec*3.28084 # convert meters to feet
    ax.plot(ot_vec, zeta_vec, '-', lw=2, c=mcolor)
    ax.plot(T['ocean_time'], zeta_vec[iot],'o', ms=20, c=mcolor)
    x0 = ot_vec[0]; x1 = ot_vec[-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta_vec.min())-5
    y1 = np.ceil(zeta_vec.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], mcolor, alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(ot_vec[0] + (ot_vec[-1]-ot_vec[0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p') + ' ' + dt_local.tzname() + '',
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c=mcolor,
        transform=ax.transAxes)
    # add day/night
    Srise = in_dict['Srise']
    Sset = in_dict['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    plt.rcdefaults()
        


