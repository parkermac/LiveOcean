# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:02:43 2016

@author: PM5

Module of basic utilities for plotting.  The goal is to make the code in
pfun less cumbersome to write and edit.

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun
import netCDF4 as nc

import numpy as np
from datetime import datetime, timedelta
import pytz

import ephem_functions as efun

if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # regular (remote, linux) version
    import matplotlib as mpl
    mpl.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import pinfo

def get_limits(Q):
    # set limits and ticklabels
    if Q['dom'] == 'full':
        Q['aa'] = []
        Q['xtl'] = range(-129,-121,2)
        Q['ytl'] = range(44,52,2)
    elif Q['dom'] == 'PS':
        Q['aa'] = [-124, -122, 47, 49.5]
        Q['xtl'] = [-123]
        Q['ytl'] = [48, 49]
        
def get_moor_info(Q):
    # set mooring info
    M = dict()
    if Q['dom'] == 'full':
        M['lon'] = -124.5
        M['lat'] = 47
        M['city'] = 'Westport'
        M['wscl'] = 20
    elif Q['dom'] == 'PS':
        M['lon'] = -122.433
        M['lat'] = 47.86
        M['city'] = 'Seattle'
        M['wscl'] = 10
    return M

def plot_time_series(ax, M, T):
    iot = zfun.find_nearest_ind(M['ot'], T['ocean_time'])
    zeta = M['zeta']*3.28084 # convert meters to feet
    ax.plot(M['ot'], zeta, '-k', lw=2)
    ax.plot(T['ocean_time'], zeta[iot],'ok', ms=20)
    x0 = M['ot'][0]; x1 = M['ot'][-1]
    x11 = T['ocean_time']
    y0 = np.floor(zeta.min())-5
    y1 = np.ceil(zeta.max())+5
    ax.fill([x0, x11, x11, x0], [y0, y0, y1, y1], 'k', alpha=.2)
    ax.axhline(c='k')
    ax.axhline(y=3, c='gray',alpha=.5)
    ax.axhline(y=-3, c='gray',alpha=.5)
    ax.text(M['ot'][0] + (M['ot'][-1]-M['ot'][0])/8, 3,
        '+3', va='center', ha='center', c='gray', weight='bold')
    ax.text(M['ot'][0] + (M['ot'][-1]-M['ot'][0])/8, -3,
        '-3', va='center', ha='center', c='gray', weight='bold')
    ax.set_yticks([])    
    ax.set_xticks([])
    dt_local = get_dt_local(T['tm'])
    ax.text(.05,.07, datetime.strftime(dt_local,'%m/%d/%Y - %I%p')+' '+dt_local.tzname(),
        transform=ax.transAxes)
    ax.text(.95,.97, 'Sea Surface Height [ft]', ha='right', va='top', style='italic', c='k',
        transform=ax.transAxes)
    # add day/night
    Srise = M['Srise']
    Sset = M['Sset']
    Srise = [dt.replace(tzinfo=None) for dt in Srise]
    Sset = [dt.replace(tzinfo=None) for dt in Sset]
    NT = len(Srise)
    for ii in range(NT):
        srise = Lfun.datetime_to_modtime(Srise[ii])
        sset = Lfun.datetime_to_modtime(Sset[ii])
        ax.plot([srise, sset],[0, 0],'-', c='orange', lw=5, alpha=.7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

def get_moor(ds0, ds1, Ldir, Q, M):
    """
    Gets a simple mooring record somewhere in the domain, to use in the plotting.
    It packs all its data in the dictionary Q, which - because of the odd
    scope behaviour of python dictionaries, is modified even as far as the
    calling function is concerned, and so does not have to be returned.
    """
    if ds0 == ds1:
        # should work for forecast or snapshot
        m_fn_list = Lfun.get_fn_list('allhours', Ldir, ds0, ds1)
    else:
        # and this is for longer time spans
        m_fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
    ot_list = []
    dt_list = []
    zeta_list = []
    uwind_list = []
    vwind_list = []
    G = zrfun.get_basic_info(m_fn_list[0], only_G=True)
    
    mi = zfun.find_nearest_ind(G['lon_rho'][0,:], M['lon'])
    mj = zfun.find_nearest_ind(G['lat_rho'][:,0], M['lat'])
    for fn in m_fn_list:
        T = zrfun.get_basic_info(fn, only_T=True)
        dt_list.append(T['tm'])
        ds = nc.Dataset(fn)
        ot_list.append(ds['ocean_time'][0])
        zeta_list.append(ds['zeta'][0,mj,mi])
        uwind_list.append(ds['Uwind'][0,mj,mi])
        vwind_list.append(ds['Vwind'][0,mj,mi])
        ds.close()
    ot = zfun.fillit(np.array(ot_list))
    zeta = zfun.fillit(np.array(zeta_list))
    uwind = zfun.fillit(np.array(uwind_list))
    vwind = zfun.fillit(np.array(vwind_list))
    uwind = zfun.filt_hanning(uwind, n=5, nanpad=False)
    vwind = zfun.filt_hanning(vwind, n=5, nanpad=False)
    M['ot'] = ot
    M['zeta'] = zeta
    M['uwind'] = uwind
    M['vwind'] = vwind
    
    # get sunrise and sunset
    sdt0 = dt_list[0] - timedelta(days=1)
    sdt1 = dt_list[-1] + timedelta(days=1)
    Srise, Sset = efun.get_sunrise_sunset(sdt0, sdt1, city=M['city'])
    # these are lists of datetimes in the UTC time zone
    M['Srise'] = Srise
    M['Sset'] = Sset

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

def add_coast(ax, dir0=Ldir['data'], color='k'):
    fn = dir0 + 'coast/coast_pnw.p'
    C = pd.read_pickle(fn)
    ax.plot(C['lon'].values, C['lat'].values, '-', color=color, linewidth=0.5)

def mask_edges(ds, fld, Q):
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    aa = Q['aa']
    clat = np.cos((np.pi/180)*(aa[2]+aa[3])/2)
    pad = .3
    fld = np.ma.masked_where(xr<aa[0]+pad, fld)
    #fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2]+pad*clat, fld)
    fld = np.ma.masked_where(yr>aa[3]-pad*clat, fld)
    return fld

def get_vlims(ds, fld, Q):
    # mask to help with choosing color limits
    xr = ds['lon_rho'][1:-1,1:-1]
    yr = ds['lat_rho'][1:-1,1:-1]
    aa = Q['aa']
    fld = np.ma.masked_where(xr<aa[0], fld)
    fld = np.ma.masked_where(xr>aa[1], fld)
    fld = np.ma.masked_where(yr<aa[2], fld)
    fld = np.ma.masked_where(yr>aa[3], fld)
    (vmin, vmax) = auto_lims(fld, vlims_fac=pinfo.vlims_fac_dict[Q['vn']])
    Q['vmin'] = vmin
    Q['vmax'] = vmax

def auto_lims(fld, vlims_fac=3):
    """
    A convenience function for automatically setting color limits.
    Input: a numpy array (masked OK)
    Output: tuple of good-guess colorscale limits for a pcolormesh plot.    
    """
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    flo = np.nanmax([np.nanmean(fld) - vlims_fac*np.nanstd(fld), np.nanmin(fld)])
    fhi = np.nanmin([np.nanmean(fld) + vlims_fac*np.nanstd(fld), np.nanmax(fld)])
    return (flo, fhi)

def get_units(ds, vn):
    try:
        units = ds[vn].units
    except AttributeError:
        units = ''
    return units

def add_bathy_contours(ax, ds, depth_levs = [200], txt=False):
    # this should work with ds being a history file Dataset, or the G dict.
    h = ds['h'][:]
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    c = 'k'
    ax.contour(lon, lat, h, depth_levs, colors=c, linewidths=0.5)
    if txt==True:
        ii = 0
        for lev in depth_levs:
            ax.text(.95, .95 - ii*.03, str(lev)+' m', c=c,
                    ha='right', transform=ax.transAxes)
            ii += 1

def add_velocity_vectors(ax, ds, fn, v_scl=3, v_leglen=0.5, nngrid=80, zlev='top', center=(.8,.05)):
    # v_scl: scale velocity vector (smaller to get longer arrows)
    # v_leglen: m/s for velocity vector legend
    xc = center[0]
    yc = center[1]
    # GET DATA
    G = zrfun.get_basic_info(fn, only_G=True)
    if zlev == 'top':
        u = ds['u'][0, -1, :, :].squeeze()
        v = ds['v'][0, -1, :, :].squeeze()
    elif zlev == 'bot':
        u = ds['u'][0, 0, :, :].squeeze()
        v = ds['v'][0, 0, :, :].squeeze()
    else:
        zfull_u = get_zfull(ds, fn, 'u')
        zfull_v = get_zfull(ds, fn, 'v')
        u = get_laym(ds, zfull_u, ds['mask_u'][:], 'u', zlev).squeeze()
        v = get_laym(ds, zfull_v, ds['mask_v'][:], 'v', zlev).squeeze()
    # ADD VELOCITY VECTORS
    # set masked values to 0
    ud = u.data; ud[u.mask]=0
    vd = v.data; vd[v.mask]=0
    # create interpolant
    import scipy.interpolate as intp
    ui = intp.interp2d(G['lon_u'][0, :], G['lat_u'][:, 0], ud)
    vi = intp.interp2d(G['lon_v'][0, :], G['lat_v'][:, 0], vd)
    # create regular grid
    aaa = ax.axis()
    daax = aaa[1] - aaa[0]
    daay = aaa[3] - aaa[2]
    axrat = np.cos(np.deg2rad(aaa[2])) * daax / daay
    x = np.linspace(aaa[0], aaa[1], int(round(nngrid * axrat)))
    y = np.linspace(aaa[2], aaa[3], int(nngrid))
    xx, yy = np.meshgrid(x, y)
    # interpolate to regular grid
    uu = ui(x, y)
    vv = vi(x, y)
    mask = uu != 0
    # plot velocity vectors
    ax.quiver(xx[mask], yy[mask], uu[mask], vv[mask],
        units='y', scale=v_scl, scale_units='y', color='k')
    ax.quiver([xc, xc] , [yc, yc], [v_leglen, v_leglen],
              [v_leglen, v_leglen],
        units='y', scale=v_scl, scale_units='y', color='k',
        transform=ax.transAxes)
    ax.text(xc+.05, yc, str(v_leglen) + ' m/s',
        horizontalalignment='left', transform=ax.transAxes)
    # note: I could also use plt.quiverkey() 
        
def add_wind(ax, M, T):
    """Add a windspeed vector with circles for scale."""
    # scl is windspeed [knots] for a 1 inch circle or arrow
    # this makes a circle 1 inch (72 points) in radius
    # and a vector 1 inch long for a windspeed of "scl" knots
    iot = zfun.find_nearest_ind(M['ot'], T['ocean_time'])
    uwind = M['uwind'][iot]
    vwind = M['vwind'][iot]
    ax.plot(M['lon'],M['lat'],'o', ms=144, mfc='None', mec='k', mew=1.5, alpha=.6)
    ax.quiver(M['lon'],M['lat'], uwind*1.94384, vwind*1.94384,
            scale=M['wscl'], scale_units='inches',
            headwidth=5,headlength=5, color='k')
            
def add_wind_text(ax, aa, M, fs):
    x0 = M['lon']
    y0 = M['lat']
    dx = aa[3] - aa[2]
    xt = x0 + dx/9
    yt = y0 - dx/30
    ax.text(xt, yt, str(M['wscl'])+' knot\nwind', c='k',
        ha='center', va='center', style='italic', size=.7*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=.3))

def add_info(ax, fn, fs=12, loc='lower_right'):
    # put info on plot
    T = zrfun.get_basic_info(fn, only_T=True)
    dt_local = get_dt_local(T['tm'])
    if loc == 'lower_right':
        ax.text(.95, .075, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .065, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
    elif loc == 'upper_right':
        ax.text(.95, .935, dt_local.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .925, dt_local.strftime('%H:%M') + ' ' + dt_local.tzname(),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)

def get_dt_local(dt, tzl='US/Pacific'):
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

def get_aa(ds):
    x = ds['lon_psi'][0,:]
    y = ds['lat_psi'][:,0]
    aa = [x[0], x[-1], y[0], y[-1]]
    return aa

def draw_box(ax, aa, linestyle='-', color='k', alpha=1, linewidth=.5, inset=0):
    aa = [aa[0]+inset, aa[1]-inset, aa[2]+inset, aa[3]-inset]
    ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]], [aa[2], aa[2], aa[3], aa[3], aa[2]],
        linestyle=linestyle, color=color, alpha=alpha, linewidth=linewidth)