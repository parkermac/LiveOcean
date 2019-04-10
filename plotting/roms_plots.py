"""
Module of plotting functions.

Each function creates, and optionally saves, a plot of fields
from a ROMS history file.

INPUT: in_dict: a tuple with information to pass to the plot, such as:
- fn: text string with the full path name of the history file to plot
- fn_out: text string with full path of output file name
- auto_vlims: a boolean governing how color limits are set
- testing: a boolean for testing (e.g. shorter, faster particle tracking)
OUTPUT: either a screen image or a graphics file

"""

# imports

# The calling function, pan_plot.py, has already put alpha on the
# path so it is on the path here too.
import Lfun
Ldir = Lfun.Lstart()
if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # fjord/boiler version
    import matplotlib as mpl
    mpl.use('Agg')
import zfun
import zrfun
from importlib import reload
import pfun; reload(pfun)
import pinfo; reload(pinfo)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta
import pandas as pd

def P_basic(in_dict):

    # START
    fig = plt.figure(figsize=(12,8)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
            #pfun.add_velocity_streams(ax, ds, in_dict['fn'])
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_ths(in_dict):
    # Plot property-property plots, like theta vs. s

    # START
    fig = plt.figure(figsize=(12,8)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    
    # make a potential density field
    import seawater as sw
    s0 = 27; s1 = 35
    th0 = 4; th1 = 20
    SS, TH = np.meshgrid(np.linspace(s0, s1, 50), np.linspace(th0, th1, 50))
    SIG = sw.dens0(SS, TH) - 1000
    
    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    h = ds['h'][:]
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    
    s = ds['salt'][:].squeeze()
    th = ds['temp'][:].squeeze()
    
    ax = fig.add_subplot(111)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Theta (deg C)')
    ax.contour(SS, TH, SIG, 20)
    
    nsub = 500
    alpha = .1
    
    mask = z > -10
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.r', alpha=alpha)
    
    mask = (z < -10) & (z > -200)
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.g', alpha=alpha)
    
    mask = z < -200
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.b', alpha=alpha)
    
    ax.set_xlim(s0, s1)
    ax.set_ylim(th0, th1)
    
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_splash(in_dict):
    # pretty picture for use in Google Earth

    # START
    fig = plt.figure(figsize=(8,12)) # (16,12) or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(111)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    
    fig.tight_layout()
    
    
    
    # FINISH
    ds.close()
    plt.show()
    plt.savefig('/Users/pm7/Desktop/splash.png', transparent=True)
        
def P_color(in_dict):
    # Code to explore cool color maps

    # START
    fig = plt.figure(figsize=(12,8)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])
    
    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = [-123.5, -122, 47,48.5]
    ii = 1
    for vn in vn_list:
        
        # things about color limits
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            vlims_fac=1
        elif vn == 'temp':
            vlims_fac=2.5
        
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                aa=aa, vlims_fac=vlims_fac)
                
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc=6) 
        fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
        elif ii == 2:
            ax.set_yticklabels('')
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_3day(in_dict):
    # Makes a plot like the Chesapeake Hypoxia Forecast Page

    # START
    fig = plt.figure(figsize=(14,6))
    fn0 = in_dict['fn']
    ds0 = nc.Dataset(fn0)
    x = ds0['lon_psi'][:]
    y = ds0['lat_psi'][:]
    
    # Get the name of the history file 3 days later
    fns = fn0.split('/')
    indir = '/'.join(fns[:-2]) + '/'
    f_string = fns[-2]
    if in_dict['list_type'] == 'snapshot':
        dstr0 = f_string[1:]
        dt0 = datetime.strptime(dstr0, '%Y.%m.%d')
        dt1 = dt0 + timedelta(days=2)
        dstr1 = datetime.strftime(dt1, '%Y.%m.%d')
        fn1 = indir + 'f' + dstr1 + '/ocean_his_0025.nc'
    elif in_dict['list_type'] == 'forecast':
        fn1 = indir + f_string + '/ocean_his_0073.nc'
    else:
        print('P_3day: unsupported list_type')
    ds1 = nc.Dataset(fn1)
    
    aa = [-124, -122, 47, 49.5]
    fs=16
    
    # PLOT CODE
    ii = 0
    for vn in range(3):
        ax = fig.add_subplot(1, 3, ii+1)
        
        vn = 'oxygen'
        nlay = 0
        laystr = 'Bottom'
        cmap=pinfo.cmap_dict[vn]
        fac=pinfo.fac_dict[vn]
        
        # don't allow auto scaling of colors
        vlims = pinfo.vlims_dict[vn]
        
        if ii == 0:
            fld0 = ds0[vn][0,nlay,1:-1,1:-1].squeeze()
            cs = ax.pcolormesh(x,y,fac*fld0,cmap=cmap, vmin=vlims[0],vmax=vlims[1])
        elif ii == 1:
            fld1 = ds1[vn][0,nlay,1:-1,1:-1].squeeze()
            cs = ax.pcolormesh(x,y,fac*fld1,cmap=cmap, vmin=vlims[0],vmax=vlims[1])
            ax.set_yticklabels([])
        elif ii == 2:
            fld2 = (fld1-fld0)/3
            if vn == 'oxygen':
                vmin=-.5
                vmax=.5
            cs = ax.pcolormesh(x,y,fac*fld2,cmap='bwr_r',vmin=vmin,vmax=vmax)
            ax.set_yticklabels([])
        cb = fig.colorbar(cs)
        cb.ax.tick_params(labelsize=fs-2)
        
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.tick_params(labelsize=fs-2)
        
        if ii==0:
            ax.set_title('Current: %s %s %s' %
                (laystr, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs)
        elif ii==1:
            ax.set_title('3 Days Later: %s %s %s' %
                (laystr, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs)
        elif ii==2:
            ax.set_title('Trend: %s %s per day' %
                (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=fs)
        ax.set_xlabel('Longitude', fontsize=fs)
        if ii == 0:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs-2, loc='upper_right')
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds0.close()
    ds1.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_basic_salish(in_dict):

    # START
    fig = plt.figure(figsize=(12,8)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = [-124, -122, 47, 49]
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'salt':
            vlims_fac = .75
        else:
            vlims_fac = 2.5
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                aa=aa, vlims_fac=vlims_fac)
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper right', borderpad=3) 
        cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        # cb.ax.tick_params(labelsize=fs1-2)
        #fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds, t_scl=.6,
                t_leglen=0.1, center=(.25, .3))
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
                v_scl=20, v_leglen=2, nngrid=100, center=(.1, .1))
            ax.set_yticklabels([])
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_basic_colvos(in_dict):

    # START
    fig = plt.figure(figsize=(12,8)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = [-122.6, -122.3, 47.35, 47.65]
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'salt':
            vlims_fac = .75
        else:
            vlims_fac = 2.5
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                aa=aa, vlims_fac=vlims_fac)
        # Inset colorbar
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(ax, width="4%", height="40%", loc='upper right', borderpad=3) 
        cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
        # cb.ax.tick_params(labelsize=fs1-2)
        #fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            # pfun.add_windstress_flower(ax, ds, t_scl=.6,
            #     t_leglen=0.1, center=(.25, .3))
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'],
                v_scl=20, v_leglen=.5, nngrid=100, center=(.1, .1))
            ax.set_yticklabels([])
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    
def P_debug(in_dict):
    # Focused on debugging

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
            # add debugging information to the plot
            # gathering info
            u = ds['u'][0,-1,:,:].squeeze()
            umax, ujmax, uimax, umin, ujmin, uimin = pfun.maxmin(u)
            v = ds['v'][0,-1,:,:].squeeze()
            vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
            eta = ds['zeta'][0,:,:].squeeze()
            emax, ejmax, eimax, emin, ejmin, eimin = pfun.maxmin(eta)
            #
            G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
            def add_info(G, ax, name, grd, vval, vj, vi, ypos, clr):
                ax.text(.98, ypos,'%s = %5.1f' % (name, vval), fontweight='bold',
                    horizontalalignment='right', transform=ax.transAxes, color=clr)
                ax.plot(G['lon_'+grd][vj,vi],G['lat_'+grd][vj,vi],'*', color=clr,
                    markeredgecolor='w',markersize=18,)
            add_info(G, ax, 'umax', 'u', umax, ujmax, uimax, .36, 'r')
            add_info(G, ax, 'umin', 'u', umin, ujmin, uimin, .33, 'orange')
            add_info(G, ax, 'vmax', 'v', vmax, vjmax, vimax, .3, 'b')
            add_info(G, ax, 'vmin', 'v', vmin, vjmin, vimin, .27, 'g')
            add_info(G, ax, 'emax', 'rho', emax, ejmax, eimax, .24, 'k')
            add_info(G, ax, 'emin', 'rho', emin, ejmin, eimin, .21, 'm')
        ii += 1
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_2D(in_dict):
    # For 2D fields.
    
    # START
    fig = plt.figure(figsize=(20,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['zeta', 'ubar', 'vbar']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict)
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('%s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        ii += 1
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_dive_vort(in_dict):
    # plots surface fields of divergence and vorticity.

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    # calculate divergence and vorticity
    u = ds['u'][0, -1, :, :]
    v = ds['v'][0, -1, :, :]
    G = zrfun.get_basic_info(in_dict['fn'], only_G=True)
    dive = ((np.diff(u, axis=1)/G['DX'][:, 1:-1])[1:-1, :]
            + (np.diff(v, axis = 0)/G['DY'][1:-1, :])[:, 1:-1])   
    x = G['lon_psi'] # matrix
    y = G['lat_psi'] # matrix
    dxp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DX'])
    dyp = zfun.interp2(x, y, G['lon_rho'], G['lat_rho'], G['DY'])
    vort = np.diff(v, axis=1)/dxp - np.diff(u, axis=0)/dyp
    aa = pfun.get_aa(ds)
    scl = 1e-4
    # panel 1
    ax = fig.add_subplot(121)
    cs = plt.pcolormesh(G['lon_psi'], G['lat_psi'], dive/scl, cmap='bwr',
                        vmin=-1, vmax=1)
    tstr = ('Surface Divergence (%0.1e $s^{-1}$)' % (scl))
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(tstr)
    pfun.add_info(ax, in_dict['fn'])
    #
    # panel 2
    ax = fig.add_subplot(122)
    cs = plt.pcolormesh(G['lon_rho'], G['lat_rho'], vort/scl, cmap='bwr',
                        vmin=-1, vmax=1)
    tstr = ('Surface Vorticity (%0.1e $s^{-1}$)' % (scl))
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title(tstr)

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    
def P_pH_Arag(in_dict):

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['PH', 'ARAG']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
        
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_willapa_oa(in_dict):

    # START
    fig = plt.figure(figsize=(12,9))
    ds = nc.Dataset(in_dict['fn'])
    
    # ** override colormaps and limits
    pinfo.vlims_dict['PH'] = (7, 8.5)
    pinfo.cmap_dict['PH'] = 'Spectral'
    pinfo.vlims_dict['ARAG'] = (0,3)
    pinfo.cmap_dict['ARAG'] = 'coolwarm_r'
    # **

    # PLOT CODE
    fs = 18
    aa = [-124.4, -123.6, 46, 47.2]
    vn_list = ['PH', 'ARAG']
    ii = 1
    for vn in vn_list:
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=0,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        cb = fig.colorbar(cs)
        cb.ax.tick_params(labelsize=fs-2)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Bottom %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
            fontsize=fs)
        ax.set_xlabel('Longitude', fontsize=fs)
        if ii == 1:
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
        elif ii == 2:
            ax.set_yticklabels('')
        ax.tick_params(labelsize=fs-2)
        ii += 1
    fig.tight_layout()
        
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_willapa_oa2(in_dict):
    # lige P_willapa_oa but for surface

    # START
    fig = plt.figure(figsize=(12,9))
    ds = nc.Dataset(in_dict['fn'])
    
    # ** override colormaps and limits
    pinfo.vlims_dict['PH'] = (7, 8.5)
    pinfo.cmap_dict['PH'] = 'Spectral'
    pinfo.vlims_dict['ARAG'] = (0,3)
    pinfo.cmap_dict['ARAG'] = 'coolwarm_r'
    # **

    # PLOT CODE
    fs = 18
    aa = [-124.4, -123.6, 46, 47.2]
    vn_list = ['ARAG', 'ARAG']
    ii = 1
    slev = 0 # first panel is the bottom
    for vn in vn_list:
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude', fontsize=fs)
        if ii == 1:
            pfun.add_bathy_contours(ax, ds, depth_levs=[10, 20, 30], txt=True)
            ax.set_title('Bottom %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
                fontsize=fs)
            ax.set_ylabel('Latitude', fontsize=fs)
            pfun.add_info(ax, in_dict['fn'], fs=fs)
        elif ii == 2:
            cb = fig.colorbar(cs)
            cb.ax.tick_params(labelsize=fs-2)
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
                fontsize=fs)
            ax.set_yticklabels('')
            ax.text(-123.9, 46.85, 'Grays\nHarbor', fontsize=fs-2, style='italic')
            ax.text(-123.86, 46.58, 'Willapa\nBay', fontsize=fs-2, style='italic')
            ax.text(-123.9, 46.04, 'Columbia\nRiver', fontsize=fs-2, style='italic')
        ax.tick_params(labelsize=fs-2)
        ii += 1
        slev = -1 # second panel is the top
    fig.tight_layout()
        
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_Carbon(in_dict):

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['alkalinity', 'TIC']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
        
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_bio(in_dict):
    
    # START
    fig = plt.figure(figsize=(20, 8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = vn_list = ['NO3', 'phytoplankton', 'alkalinity', 'TIC']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn in []:
            slev = 0
            ttag = 'Bottom'
        else:
            slev = -1
            ttag = 'Surface'
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        if True:
            ax.axis(pfun.get_aa(ds))
        else:
            aa = [-127.2, -123.8, 45.5, 49.8]
            ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('%s %s %s' % (ttag,pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
        elif ii == 2:
            pfun.add_windstress_flower(ax, ds, center=(.2,.25))
        ii += 1
        
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer(in_dict):

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    z_level = -30
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        # cb.formatter.set_useOffset(False)
        # cb.update_ticks()
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=1, v_leglen=0.5,
                                      nngrid=80, zlev=z_level)
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer_JdF_Eddy(in_dict):

    # START
    fig = plt.figure(figsize=(24,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    aa = [-125.8, -124.2, 48, 49]
    pinfo.vlims_dict['salt'] = (32,34)
    pinfo.vlims_dict['temp'] = (6,10)
    
    vn_list = ['salt', 'temp']
    z_level = -50
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        # if in_dict['auto_vlims']:
        #     pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        # cb.formatter.set_useOffset(False)
        # cb.update_ticks()
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_bathy_contours(ax, ds, depth_levs = [50, 100, 200], txt=True)
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds, t_scl=0.6, t_leglen=0.1, center=(.15,.75))
        if ii == 2:
            # pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=5, v_leglen=0.2,
            #                           nngrid=30, zlev=z_level)
            pass
        ii += 1

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect(in_dict):
    # plots a map and a section (distance, z)
    
    # START
    fig = plt.figure(figsize=(24,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    vn = 'phytoplankton'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -3500
        #x = np.linspace(lon.min(), -124, 500)
        if True:
            x = np.linspace(lon.min(), lon.max(), 500)
            y = 47 * np.ones(x.shape)
        else:
            y = np.linspace(lat.min(), lat.max(), 500)
            x = -126 * np.ones(y.shape)
    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] + 'tracks_new/'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path + track
            # get the track to interpolate onto
            pdict = pickle.load(open(track_fn, 'rb'))
            xx = np.concatenate((xx,pdict['lon_poly']))
            yy = np.concatenate((yy,pdict['lat_poly']))
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect2(in_dict):
    # Plots a map and a section (distance, z) for 2 variables.
    
    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    #tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
    tracks = ['Line_jdf_v0.p', 'Line_HC_thalweg_long.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))

    # PLOTTING
    #vn_list = ['NO3', 'phytoplankton']
    #vn_list = ['NO3', 'oxygen']
    vn_list = ['PH', 'ARAG']
    # override colors
    pinfo.vlims_dict['NO3'] = (0,40)
    pinfo.vlims_dict['phytoplankton'] = (0,20)
    pinfo.vlims_dict['oxygen'] = (0,8)
    pinfo.vlims_dict['PH'] = (7, 8)
    pinfo.vlims_dict['ARAG'] = (0, 3)
    counter = 0
    for vn in vn_list:
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

        # map with section line
        if counter == 0:
            ax = fig.add_subplot(2, 3, 1)
        elif counter == 1:
            ax = fig.add_subplot(2, 3, 4)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        pfun.add_coast(ax)
        ax.axis([-126, -122, 47, 50])
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_ylabel('Latitude')
        if counter == 1:
            ax.set_xlabel('Longitude')
        # add section track
        ax.plot(x, y, '-w', linewidth=2)
        ax.plot(x, y, '-k', linewidth=.5)
        ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
            markeredgecolor='r', markeredgewidth=2)
        #
        # section
        if counter == 0:
            ax = fig.add_subplot(2, 3, (2, 3))
        elif counter == 1:
            ax = fig.add_subplot(2, 3, (5, 6))
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(dist.min(), dist.max())
        ax.set_ylim(zdeep, 5)
        sf = pinfo.fac_dict[vn] * v3['sectvarf']
        # set section color limits
        vmin = pinfo.vlims_dict[vn][0]
        vmax = pinfo.vlims_dict[vn][1]
        # plot section
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                           vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
        fig.colorbar(cs)
        cs = ax.contour(v3['distf'], v3['zrf'], sf,
            np.linspace(np.floor(vmin), np.ceil(vmax), 20),
            colors='k', linewidths=0.5)
        if counter == 1:
            ax.set_xlabel('Distance (km)')
            pfun.add_info(ax, in_dict['fn'])
        ax.set_ylabel('Z (m)')
        ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        
        counter += 1
        
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot(in_dict):
    # Plot salinity maps and section, with forcing time-series.
    # Super clean design.

    vn = 'salt'
    vlims = (28.5, 33) # full map
    vlims2 = (22, 31) # PS map
    vlims3 = (29, 32) # PS section
    cmap = 'jet'

    # get model fields
    fn = in_dict['fn']
    ds = nc.Dataset(fn)

    # get forcing fields
    ffn = Ldir['LOo'] + 'superplot/forcing_cas4_v2_lo6biom_2017.p'
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] + 'tracks_new/'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'][:]
    lat = ds['lat_psi'][:]
    v =ds[vn][0, -1, 1:-1, 1:-1]
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nSalinity\n'
        + datetime.strftime(T['tm'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    
    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['tm'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(.4, 1.7)
    ax.set_axis_off()
    
    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)    
    ax.set_ylim(-5,20)
    ax.set_axis_off()
    
    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()    

def P_superplot2(in_dict):
    # Plot phytoplankton maps and section, with forcing time-series.
    # Super clean design.  Simpler than P_superplot.

    vn = 'phytoplankton'
    vlims = (0, 20)
    cmap = 'viridis'
    sect_color = 'orange'
    up_color = 'gray'
    down_color = 'gray'
    now_color = 'brown'
    fs = 16 # fontsize
    aa = [-123.3, -122.1, 47.01, 48.4]

    # get model fields
    fn = in_dict['fn']
    ds = nc.Dataset(fn)

    # get forcing fields
    ffn = Ldir['LOo'] + 'superplot/forcing_cas4_v2_lo6biom_2017.p'
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] + 'tracks_new/'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path + track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'][:]
    lat = ds['lat_psi'][:]
    v =ds[vn][0, -1, 1:-1, 1:-1]
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    pfun.draw_box(ax, aa, color=sect_color, alpha=1, linewidth=5, inset=.01)
    # labels
    ax.text(.03, .02, 'LiveOcean\nPhytoplankton\n'
        + datetime.strftime(T['tm'], '%Y'), fontsize=fs, color='w',
        transform=ax.transAxes, horizontalalignment='left',
        fontweight='bold')

    # PS map
    ax =  plt.subplot2grid((3,3), (1,2), rowspan=2)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1],
        cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color=sect_color, alpha=1, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle=':', color=sect_color, linewidth=3)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')

    # Section
    ax =  plt.subplot2grid((3,3), (0,1), colspan=2)
    ax.plot(dist, v2['zeta']+5, linestyle=':', color=sect_color, linewidth=3)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 30)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=0, vmax=20, cmap=cmap)
    # labels
    ax.text(1, .2, 'Puget Sound Section', fontsize=fs, color=sect_color,
        transform=ax.transAxes, horizontalalignment='right', fontweight='bold')
    ax.set_axis_off()

    # get the day
    tm = T['tm'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    
    # Wind
    alpha=1
    ax =  plt.subplot2grid((3,3), (1,1), rowspan=2)
    #ax = fig.add_subplot(336)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color=down_color, alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color=up_color, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color=now_color, markersize=14)
    # labels
    ax.text(0, .85, 'Wind from South:\nDecreases Nutrients at Coast', transform=ax.transAxes,
        color=down_color, alpha=alpha, fontsize=fs, fontweight='bold')
    ax.text(0, .25, 'Wind from North:\nBrings Nutrients to Surface at Coast', transform=ax.transAxes,
        color=up_color, alpha=alpha, fontsize=fs, fontweight='bold', verticalalignment='bottom')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.2, .25)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = now_color
    if tm.month in [4, 5, 6]:
        clist[1] = now_color
    if tm.month in [7, 8, 9]:
        clist[2] = now_color
    if tm.month in [10, 11, 12]:
        clist[3] = now_color
    ypos = .1
    tfs = fs + 4
    ax.text(0, ypos, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=tfs, horizontalalignment='left', fontweight='bold')
    ax.text(.39, ypos, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=tfs, horizontalalignment='center', fontweight='bold')
    ax.text(.67, ypos, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=tfs, horizontalalignment='center', fontweight='bold')
    ax.text(1, ypos, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=tfs, horizontalalignment='right', fontweight='bold')
    
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_hc(in_dict):
    # plots a map and a section (distance, z)
    
    # START
    fig = plt.figure(figsize=(24,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if vn == 'salt':
        vlims_fac = 1
    else:
        vlims_fac = 2.5
    
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    zdeep = -120
    track = 'HC_north.p'
    track_fn = tracks_path + track
    # get the track to interpolate onto
    pdict = pickle.load(open(track_fn, 'rb'))
    xx = pdict['lon_poly']
    yy = pdict['lat_poly']
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    aa = [-123.2, -122.4, 47.2, 48.2]
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], aa=aa, vlims_fac=vlims_fac)
    fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    #pfun.add_windstress_flower(ax, ds)
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf, vlims_fac=vlims_fac)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_willapa(in_dict):
    # Focus on Willapa
    
    # override
    in_dict['auto_vlims'] = False
    
    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    vn = 'phytoplankton'
    #vn = 'salt'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    else:
        pinfo.vlims_dict['sect_'+vn] = pinfo.vlims_dict[vn]
        
    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    track = 'Line_willapa1.p'
    zdeep = -30
    xx = np.array([])
    yy = np.array([])
    track_fn = tracks_path + track
    # get the track to interpolate onto
    pdict = pickle.load(open(track_fn, 'rb'))
    xx = np.concatenate((xx,pdict['lon_poly']))
    yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    
    # full map
    ax = fig.add_subplot(1, 4, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
        
    # focus map
    ax = fig.add_subplot(1, 4, 2)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    ax.axis([-124.4, -123.6, 46, 47.2])
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 4, (3, 4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect_columbia(in_dict):
    # Focus on Columbia River
    
    # override
    in_dict['auto_vlims'] = False
    
    # START
    fig = plt.figure(figsize=(18,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    #vn = 'NO3'
    #vn = 'phytoplankton'
    vn = 'salt'
    # we allow for the possibility of using different color scales
    # on the map and section for the same varible, and these follow
    # the general scheme that the default is for them to be chosen
    # automatically based on the first plot in a series.
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
        pinfo.vlims_dict['sect_'+vn] = ()
    else:
        pinfo.vlims_dict[vn] = (0,35)
        pinfo.vlims_dict['sect_'+vn] = (0,35)
        
    #
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    tracks_path = Ldir['data'] + 'tracks_new/'
    track = 'CR_thalweg.p'
    zdeep = -30
    xx = np.array([])
    yy = np.array([])
    track_fn = tracks_path + track
    # get the track to interpolate onto
    pdict = pickle.load(open(track_fn, 'rb'))
    xx = np.concatenate((xx,pdict['lon_poly']))
    yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    
    # full map
    ax = fig.add_subplot(1, 4, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    pfun.add_info(ax, in_dict['fn'])
    pfun.add_windstress_flower(ax, ds)
    
    aa = [-124.5, -123.3, 45.7, 47.2]
    ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]],
        [aa[2], aa[2], aa[3], aa[3], aa[2]], '-c', linewidth=2)
        
    # focus map
    ax = fig.add_subplot(1, 4, 2)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    #
    # section
    ax = fig.add_subplot(1, 4, (3, 4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # set section color limits
    svlims = pinfo.vlims_dict['sect_'+vn]
    if len(svlims) == 0:
        svlims = pfun.auto_lims(sf)
        pinfo.vlims_dict['sect_'+vn] = svlims
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], sf,
        np.linspace(np.floor(svlims[0]), np.ceil(svlims[1]), 20),
        colors='k', linewidths=0.5)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    
def P_sectA(in_dict):
    # plots a map and several sections
    # designed for analytical runs, like aestus1
    # used for MacCready et al. (2018 JPO) variance paper

    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])
    
    # PLOTTING
    vn = 'salt'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    # map and section lines
    ax1 = fig.add_subplot(3,1,1)
    cs = pfun.add_map_field(ax1, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn])
    vlims = pinfo.vlims_dict[vn]
    fig.colorbar(cs)
    ax1.axis([-.5, 1, 44.8, 45.2])
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    pfun.add_info(ax1, in_dict['fn'], fs=9)
    #
    # thalweg section
    x = np.linspace(-0.5, 1, 500)
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    ax = fig.add_subplot(3,1,2)
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-25, 2)
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                       vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.linspace(1, 35, 35), colors='k', linewidths=.5,)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(vn)
    # add line to map plot
    ax1.plot(x, y, '-k', linewidth=2)
    ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
        markeredgecolor='k', markeredgewidth=2)
    # channel cross-sections
    for ii in range(3):
        y = np.linspace(44.9, 45.1, 200)
        x = 0.3*ii * np.ones(y.shape)
        v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
        ax = fig.add_subplot(3,3,ii+7)
        ax.plot(dist, v2['zbot'], '-k', linewidth=2)
        ax.plot(dist, v2['zeta'], '-b', linewidth=1)
        ax.set_xlim(0, 22.5)
        ax.set_ylim(-25, 2)
        cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
            np.linspace(1, 35, 35), colors='k', linewidths=.5,)
        ax.set_xlabel('Distance (km)')
        if ii==0:
            ax.set_ylabel('Z (m)')
        # add line to map plot
        ax1.plot(x, y, '-k', linewidth=2)
        ax1.plot(x[idist0], y[idist0], 'ok', markersize=10, markerfacecolor='w',
            markeredgecolor='k', markeredgewidth=2)
    #
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tracks_MERHAB(in_dict):
    # Use trackfun_1 to create surface drifter tracks for MERHAB.
    # It automatically makes tracks for as long as there are
    # hours in the folder.
    
    # START
    fig = plt.figure(figsize=(12, 8))
    ds = nc.Dataset(in_dict['fn'])

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle

    # need to get Ldir so that we can get a full list
    # of the history files for this day,
    # which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    # and estimate the number of days,
    ndays = round(len(fn_list_full)/24)
    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])

    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    T0 = zrfun.get_basic_info(fn_list[0], only_T=True)

    if len(fn_list) == 2:
        # only do the tracking at the start
    
        # standard MERHAB version
        nyp = 7
        x0 = -126
        x1 = -125
        y0 = 48
        y1 = 49
        mlr = np.pi*(np.mean([y0, y1]))/180
        xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
        lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        latvec = np.linspace(y0, y1, nyp)
        lonmat_1, latmat_1 = np.meshgrid(lonvec, latvec)
        x0 = -125
        x1 = -124
        y0 = 44
        y1 = 45
        mlr = np.pi*(np.mean([y0, y1]))/180
        xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
        lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        latvec = np.linspace(y0, y1, nyp)
        lonmat_2, latmat_2 = np.meshgrid(lonvec, latvec)
        lonmat = np.concatenate((lonmat_1.flatten(), lonmat_2.flatten()))
        latmat = np.concatenate((latmat_1.flatten(), latmat_2.flatten()))
        #
        plon00 = lonmat.flatten()
        plat00 = latmat.flatten()
        pcs00 = np.array([-.05]) # unimportant when surface=True
        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        # DO THE TRACKING
        import time
        tt0 = time.time()
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': 1, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        # load the output on all subsequent days
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE

    # panel 1
    ax = fig.add_subplot(121)
    vn = 'salt'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            alpha=0.5)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    fs1 = 14
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    tstr = 'Surface ' + pinfo.tstr_dict[vn] +' and ' + str(ndays) + ' day Tracks'
    ax.set_title(tstr, fontsize=fs1)
    # add info
    fs = fs1 - 6
    ax.text(.98, .10, T0['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right' , verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .07, 'to ' + T['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .04, 'UTC',
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)
    # add the tracks
    c_start = 'w'; s_start = 4
    c_end = 'r'; s_end = 6
    ntt = len(fn_list) -1
    ax.plot(P['lon'][:ntt,:], P['lat'][:ntt,:], '-k', linewidth=2)
    ax.plot(P['lon'][0,:],P['lat'][0,:],'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k')
    ax.plot(P['lon'][ntt,:],P['lat'][ntt,:],'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k')
    # add info about the tracks
    x0 = .8; x1 = .9
    y0 = .25; y1 = .3
    ax.plot([x0, x1], [y0, y1], '-k', linewidth=2, transform=ax.transAxes)
    ax.plot(x0, y0,'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    ax.plot(x1, y1,'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    # and some labels
    ax.text(x0+.02, y0, 'start', horizontalalignment='left',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
    ax.text(x1-.02, y1, 'end', horizontalalignment='right',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
    #
    ax = fig.add_subplot(122)
    vn = 'phytoplankton'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    tstr = 'Surface ' + pinfo.tstr_dict[vn]
    if gtx_list[0]=='cas4':
        mask_sal=False
    else:
        mask_sal=True
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            do_mask_salish=mask_sal)
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=False)
    pfun.add_coast(ax)
    pfun.add_windstress_flower(ax, ds)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        fontsize=fs1)
    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_merhab2(in_dict):

    # START
    fig = plt.figure(figsize=(14, 12))
    ds = nc.Dataset(in_dict['fn'])

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle

    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    # and estimate the number of days,
    ndays = round(len(fn_list_full)/24)
    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])

    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    T0 = zrfun.get_basic_info(fn_list[0], only_T=True)
    
    # Create initial positions for tracking
    nyp = 7
    x0 = -126; x1 = -125; y0 = 48; y1 = 49
    clat_1 = np.cos(np.pi*(np.mean([y0, y1]))/180)
    xyRatio = clat_1 * (x1 - x0) / (y1 - y0)
    lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
    latvec = np.linspace(y0, y1, nyp)
    lonmat_1, latmat_1 = np.meshgrid(lonvec, latvec)
    #
    x0 = -125; x1 = -124; y0 = 44; y1 = 45
    clat_2 = np.cos(np.pi*(np.mean([y0, y1]))/180)
    xyRatio = clat_2 * (x1 - x0) / (y1 - y0)
    lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
    latvec = np.linspace(y0, y1, nyp)
    lonmat_2, latmat_2 = np.meshgrid(lonvec, latvec)
    lonmat = np.concatenate((lonmat_1.flatten(), lonmat_2.flatten()))
    latmat = np.concatenate((latmat_1.flatten(), latmat_2.flatten()))
    #
    plon00 = lonmat.flatten(); plat00 = latmat.flatten()
    pcs00 = np.array([-.05])

    if (len(fn_list)==2) and (in_dict['testing']==False):
        # only do the tracking at the start
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        # DO THE TRACKING
        import time
        tt0 = time.time()
        ndiv = 1
        print('** using ndiv = 1 **')
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': ndiv, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        # load the output on all subsequent days
        if in_dict['testing'] == True:
            outdir = in_dict['outdir']
        else:
            fo = in_dict['fn_out']
            outdir = fo[:fo.rindex('/')+1]
        out_fn = outdir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE
    fs1 = 18
    
    # panel 1
    ax = fig.add_subplot(121)
    vn = 'salt'
    vlims = (28.5,33)
    pinfo.vlims_dict['salt'] = vlims
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    pfun.add_windstress_flower(ax, ds, fs=fs1-4)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    tstr = ('Surface Salinity ['+str(vlims[0])+'-'+str(vlims[1])+']')
    ax.set_title(tstr, fontsize=fs1)
    # add info
    pfun.add_info(ax, in_dict['fn'], fs=fs1)
    # add the tracks
    ntt = len(fn_list)
    ax.plot(P['lon'][:,:], P['lat'][:,:], '-w', linewidth=2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1-3)
    ax.text(.05,.5, 'Three Day Tracks', color='w', fontsize=fs1-2, transform=ax.transAxes)
    
    #
    ax1 = fig.add_subplot(222)
    pad = .7
    aa1 = [ -126 - pad/clat_1 - pad/3, -125 + pad/clat_1 - pad/3, 48 - pad, 49 + pad]
    pfun.draw_box(ax, aa1, linestyle='-', color='c', alpha=1, linewidth=3)
    pfun.draw_box(ax1, aa1, linestyle='-', color='c', alpha=1, linewidth=5, inset=.01)
    pfun.add_coast(ax1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    # grid of initial positions
    pnr, pnc = lonmat_1.shape
    for jj in range(pnc):
        ax1.plot([lonmat_1[0,jj],lonmat_1[-1,jj]],[latmat_1[0,jj],latmat_1[-1,jj]],'-c')
    for ii in range(pnr):
        ax1.plot([lonmat_1[ii,0],lonmat_1[ii,-1]],[latmat_1[ii,0],latmat_1[ii,-1]],'-c')
    # add the tracks
    for ii in range(ntt):
        if ii > ntt-24:
            mksz = (ii - ntt + 24)/2
            mkco = 'red'
        else:
            mksz = 1
            mkco = 'gray'
        ax1.plot(P['lon'][ii,:], P['lat'][ii,:], 'o', color=mkco, markersize=mksz, alpha=.3)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = fs1-3)
    ax1.text(.05,.05, 'Juan de Fuca Eddy\n24 Hour Tracks', color='r', fontsize=fs1-2, transform=ax1.transAxes)
    
    ax2 = fig.add_subplot(224)
    aa2 = [ -125 - pad/clat_2 - pad/3, -124 + pad/clat_2 - pad/3, 44 - pad, 45 + pad]
    pfun.draw_box(ax, aa2, linestyle='-', color='g', alpha=1, linewidth=3)
    pfun.draw_box(ax2, aa2, linestyle='-', color='g', alpha=1, linewidth=5, inset=.01)
    pfun.add_coast(ax2)
    ax2.axis(aa2)
    pfun.dar(ax2)    
    # grid of initial positions
    pnr, pnc = lonmat_2.shape
    for jj in range(pnc):
        ax2.plot([lonmat_2[0,jj],lonmat_2[-1,jj]],[latmat_2[0,jj],latmat_2[-1,jj]],'-g')
    for ii in range(pnr):
        ax2.plot([lonmat_2[ii,0],lonmat_2[ii,-1]],[latmat_2[ii,0],latmat_2[ii,-1]],'-g')
    # add the tracks
    for ii in range(ntt):
        if ii > ntt-24:
            mksz = (ii - ntt + 24)/2
            mkco = 'red'
        else:
            mksz = 1
            mkco = 'gray'
        ax2.plot(P['lon'][ii,:], P['lat'][ii,:], 'o', color=mkco, markersize=mksz, alpha=.3)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = fs1-3)
    ax2.text(.05,.05, 'Heceta Bank\n24 Hour Tracks', color='r', fontsize=fs1-2, transform=ax2.transAxes)
        
    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tracks_full(in_dict):
    # Use trackfun_1 to create surface drifter tracks for a full-domain plot.
    # It automatically makes tracks for as long as there are
    # hours in the folder.
    
    # START
    fig = plt.figure(figsize=(8,12))
    ds = nc.Dataset(in_dict['fn'])

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle

    # need to get Ldir so that we can get a full list
    # of the history files for this day,
    # which means unpacking gtagex
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    # and estimate the number of days,
    ndays = round(len(fn_list_full)/24)
    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])

    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    T0 = zrfun.get_basic_info(fn_list[0], only_T=True)

    if len(fn_list) == 2:
        # only do the tracking at the start
        nyp = 70
        x0 = -127.1
        x1 = -122
        y0 = 42.3
        y1 = 50
        mlr = np.pi*(np.mean([y0, y1]))/180
        xyRatio = np.cos(mlr) * (x1 - x0) / (y1 - y0)
        lonvec = np.linspace(x0, x1, (nyp * xyRatio).astype(int))
        latvec = np.linspace(y0, y1, nyp)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        # 2019.03.11 limit to shallow water
        hh = zfun.interp2(lonmat, latmat, G['lon_rho'], G['lat_rho'], G['h'])
        hmask = hh <= 400
        lonmat = lonmat[hmask]
        latmat = latmat[hmask]
        #
        plon00 = lonmat.flatten()
        plat00 = latmat.flatten()
        pcs00 = np.array([-.05]) # unimportant when surface=True
        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        
        # DO THE TRACKING
        import time
        tt0 = time.time()
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': 1, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        # load the output on all subsequent days
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE

    # panel 1
    ax = fig.add_subplot(111)
    vn = 'salt'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            alpha=0.5)
    fig.colorbar(cs)
    #pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    fs1 = 14
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    tstr = 'Surface ' + pinfo.tstr_dict[vn] +' and ' + str(ndays) + ' day Tracks'
    ax.set_title(tstr, fontsize=fs1)
    # add info
    fs = fs1 - 2
    ax.text(.98, .10, T0['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right' , verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .07, 'to ' + T['tm'].strftime('%Y-%m-%d %H:%M'),
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.98, .04, 'UTC',
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=fs)
    ax.text(.06, .04, fn.split('/')[-3],
        verticalalignment='bottom', transform=ax.transAxes,
        rotation='vertical', fontsize=fs)
    # add the tracks
    t_lw = .5
    c_start = 'w'; s_start = 1
    c_end = 'r'; s_end = 3
    ntt = len(fn_list) -1
    ax.plot(P['lon'][:ntt,:], P['lat'][:ntt,:], '-k', linewidth=t_lw)
    ax.plot(P['lon'][0,:],P['lat'][0,:],'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k')
    ax.plot(P['lon'][ntt,:],P['lat'][ntt,:],'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k')
    # add info about the tracks
    x0 = .8; x1 = .9
    y0 = .25; y1 = .3
    ax.plot([x0, x1], [y0, y1], '-k', linewidth=t_lw, transform=ax.transAxes)
    ax.plot(x0, y0,'o'+c_start,
            markersize=s_start, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    ax.plot(x1, y1,'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k', transform=ax.transAxes)
    # and some labels
    ax.text(x0+.02, y0, 'start', horizontalalignment='left',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
    ax.text(x1-.02, y1, 'end', horizontalalignment='right',
            verticalalignment='center', fontstyle='italic', transform=ax.transAxes)
    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tracks_ps(in_dict):
    # Use trackfun_1 to create surface drifter tracks for Puget Sound.
    # It automatically makes tracks for as long as there are
    # hours in the folder.
    
    # START
    fig = plt.figure(figsize=(16,12))#(10,12))

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
            
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    
    if in_dict['testing'] == True:
        fn_list_full = fn_list_full[:11]

    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])
    
    if len(fn_list) == 2:
        # only do the tracking at the start

        # Specific release locations
        rloc_dict = {'Seattle': (-122.36, 47.6),
                    'Tacoma': (-122.44, 47.285)}
        count = 0
        rnp = 1000
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]
            ry = rloc_dict[rloc][1]
            #print('%s %0.3f %0.3f' % (rloc, rx, ry))
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy()
                plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        #
        pcs00 = np.array([-.05])
        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        # DO THE TRACKING
        import time
        tt0 = time.time()
        # NOTE: we use at least ndiv=12 to get advection right in places
        # like Tacoma Narrows, but we reduce it for testing.
        if (Ldir['lo_env'] == 'pm_mac') or (in_dict['testing'] == True):
            ndiv = 1
            print('** using ndiv = 1 **')
        else:
            ndiv = 12
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': ndiv, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        # load the output on all subsequent days
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE
    fs1 = 18
    vn = 'salt'
    vlims_fac = 1
    # if in_dict['auto_vlims']:
    #     pinfo.vlims_dict[vn] = ()
    # override
    pinfo.vlims_dict['salt'] = (28.5,33)
    
    # panel 2: do this first so it controls the color limits (obsolete, right?)
    ax = fig.add_subplot(122)
    aa = [-123.5, -122, 47, 48.5]
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            alpha=1, aa=aa, vlims_fac=vlims_fac)
    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="4%", height="40%", loc=6) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    cb.ax.tick_params(labelsize=fs1-2)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    nhour = 3
    tstr = ('Puget Sound: Surface ' + pinfo.tstr_dict[vn] +
            ' and ' + str(nhour) + ' hour Particle Tracks')
    ax.set_title(tstr, fontsize=fs1)
    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1-2)
    # add the tracks
    t_lw = .5
    c_end = 'k'; s_end = 1
    ntt = len(fn_list) -1
    ntt0 = ntt-nhour
    if ntt0<0:
        ntt0 = 0
    ax.plot(P['lon'][ntt0:ntt+1,:], P['lat'][ntt0:ntt+1,:], '-k', linewidth=t_lw)
    ax.plot(P['lon'][ntt,:],P['lat'][ntt,:],'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k')
    
    # panel 1
    ax = fig.add_subplot(121)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    #fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_title('Full Domain: Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
        fontsize=fs1)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    pfun.add_info(ax, in_dict['fn'], fs=fs1)
    pfun.add_windstress_flower(ax, ds, fs=fs1-2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1-2)
        

    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_tracks_fhl(in_dict):
    # Use trackfun_1 to create surface drifter tracks for FHL.
    # It automatically makes tracks for as long as there are
    # hours in the folder.
    
    # START
    fig = plt.figure(figsize=(16,12))#(10,12))

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
            
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    
    if in_dict['testing'] == True:
        fn_list_full = fn_list_full[:11]

    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])
    
    if len(fn_list) == 2:
        # only do the tracking at the start

        # Specific release locations
        rloc_dict = {'FHL': (-123.0, 48.545)}
        count = 0
        rnp = 1000
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]
            ry = rloc_dict[rloc][1]
            #print('%s %0.3f %0.3f' % (rloc, rx, ry))
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy()
                plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        #
        pcs00 = np.array([-.05])
        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        # DO THE TRACKING
        import time
        tt0 = time.time()
        # NOTE: we use at least ndiv=12 to get advection right in places
        # like Tacoma Narrows, but we reduce it for testing.
        if (in_dict['testing'] == True):
            ndiv = 1
            print('** using ndiv = 1 **')
        else:
            ndiv = 12
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': ndiv, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
    else:
        # load the output on all subsequent days
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        P = pickle.load(open(out_fn, 'rb'))

    # PLOT CODE
    fs1 = 18
    vn = 'salt'
    vlims_fac = 1
    # if in_dict['auto_vlims']:
    #     pinfo.vlims_dict[vn] = ()
    # override
    pinfo.vlims_dict['salt'] = (23,33)
    
    # panel 2: do this first so it controls the color limits
    ax = fig.add_subplot(122)
    aa = [-123.5, -122.2, 47.8, 49.3]
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            alpha=1, aa=aa, vlims_fac=vlims_fac)
    # # Inset colorbar
    # from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # cbaxes = inset_axes(ax, width="4%", height="40%", loc=7, borderpad=1)
    # cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    # cb.ax.tick_params(labelsize=fs1-2)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    #pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=3, inset=.01)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    nhour = 12
    tstr = ('San Juan Islands: Surface ' + pinfo.tstr_dict[vn] +
            ' (' + str(pinfo.vlims_dict[vn][0]) + '-' + str(pinfo.vlims_dict[vn][1]) + ')' +
            ' and ' + str(nhour) + ' hour Tracks')
    ax.set_title(tstr, fontsize=fs1)
    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1-2)
    # add the tracks
    t_lw = .5
    c_end = 'k'; s_end = 1
    ntt = len(fn_list) -1
    ntt0 = ntt-nhour
    if ntt0<0:
        ntt0 = 0
    ax.plot(P['lon'][ntt0:ntt+1,:], P['lat'][ntt0:ntt+1,:], '-k', linewidth=t_lw)
    ax.plot(P['lon'][ntt,:],P['lat'][ntt,:],'o'+c_end,
            markersize=s_end, alpha = 1, markeredgecolor='k')
    
    # panel 1
    ax = fig.add_subplot(121)
    # override
    pinfo.vlims_dict['salt'] = (28.5,33)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
    # # Inset colorbar
    # from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # cbaxes = inset_axes(ax, width="4%", height="40%", loc=7, borderpad=1)
    # cb = fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    # cb.ax.tick_params(labelsize=fs1-2)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=3, inset=.01)
    ax.set_title('Full Domain: Surface %s %s (%s-%s)' % (pinfo.tstr_dict[vn],
        pinfo.units_dict[vn],str(pinfo.vlims_dict[vn][0]),str(pinfo.vlims_dict[vn][1])),
        fontsize=fs1)
    ax.set_xlabel('Longitude', fontsize=fs1)
    ax.set_ylabel('Latitude', fontsize=fs1)
    pfun.add_info(ax, in_dict['fn'], fs=fs1)
    pfun.add_windstress_flower(ax, ds, fs=fs1-2)
    ax.tick_params(axis = 'both', which = 'major', labelsize = fs1-2)
        

    #
    fig.tight_layout()

    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_tracks_ps_debug(in_dict):
    # Figure out why P_tracks_ps sometimes gets ot = 0 in hour 0.
    
    # START

    # TRACKS
    import os
    import sys
    pth = os.path.abspath('../tracker')
    if pth not in sys.path:
        sys.path.append(pth)
    import trackfun_1 as tf1
    import pickle
    fn = in_dict['fn']
    gtagex = fn.split('/')[-3]
    gtx_list = gtagex.split('_')
    import Lfun
    Ldir = Lfun.Lstart(gtx_list[0], gtx_list[1])
    # make a list of history files in this folder,
    in_dir = fn[:fn.rindex('/')+1]
    fn_list_raw = os.listdir(in_dir)
    fn_list = []
    for item in fn_list_raw:
        if 'ocean_his' in item:
            fn_list.append(in_dir + item)
    fn_list.sort()
    
            
    # save the full list to use with the tracking code 
    fn_list_full = fn_list.copy()
    
    if in_dict['testing'] == True:
        fn_list_full = fn_list_full[:2]

    # then trim fn_list to end at the selected hour for this plot
    fn_list = fn_list[:fn_list.index(fn)+1]
    
    # and use the CURRENT file for the map field overlay
    ds = nc.Dataset(in_dict['fn'])
    
    if len(fn_list) == 2:
        # only do the tracking at the start

        # Specific release locations
        rloc_dict = {'Seattle': (-122.36, 47.6),
                    'Tacoma': (-122.44, 47.285)}
        count = 0
        rnp = 1000
        for rloc in rloc_dict:
            rx = rloc_dict[rloc][0]
            ry = rloc_dict[rloc][1]
            #print('%s %0.3f %0.3f' % (rloc, rx, ry))
            px00 = rx + 0.01*np.random.randn(rnp)
            py00 = ry + 0.01*np.random.randn(rnp)
            if count == 0:
                plon00 = px00.copy()
                plat00 = py00.copy()
            else:
                plon00 = np.concatenate((plon00, px00.copy()))
                plat00 = np.concatenate((plat00, py00.copy()))
            count += 1
        #
        pcs00 = np.array([-.05])
        # make the full IC vectors, which will have equal length
        # (one value for each particle)
        NSP = len(pcs00)
        NXYP = len(plon00)
        plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
        plon0 = plon0.flatten()
        plat0 = plat0.flatten()
        pcs0 = pcs0.flatten()
        # DO THE TRACKING
        import time
        tt0 = time.time()
        # NOTE: we use at least ndiv=12 to get advection right in places
        # like Tacoma Narrows, but we reduce it for testing.
        if (Ldir['lo_env'] == 'pm_mac') or (in_dict['testing'] == True):
            ndiv = 1
            #print('** using ndiv = 1 **')
        else:
            ndiv = 12
        TR = {'3d': False, 'rev': False, 'turb': False,
            'ndiv': ndiv, 'windage': 0}
        P = tf1.get_tracks(fn_list_full, plon0, plat0, pcs0, TR,
                           trim_loc=True)
        #print('  took %0.1f seconds' % (time.time() - tt0))
        # and store the output
        fo = in_dict['fn_out']
        out_dir = fo[:fo.rindex('/')+1]
        out_fn = out_dir + 'tracks.p'
        pickle.dump(P, open(out_fn, 'wb'))
        
        #print('\nLooking at P')
        rot = P['ot']
        if rot[0] == 0:
            print('Error Found!')
            time_dt = Lfun.modtime_to_datetime(item)
            print('%s %s' % (str(item),str(time_dt)))
            sys.exit()
        
        # x = P['lon']
        # for item in x[:,0]:
        #     print('x = %0.4f' % (item))
        # print('')
        
    else:
        #print('len(fn_list) = %d' % (len(fn_list)))
        sys.exit()






