"""
Module of plotting functions.

Each function creates, and optionally saves, a plot of fields
from a ROMS history file.

INPUT: in_dict: a tuple with information to pass to the plot, such as:
- fn: text string with the full path name of the history file to plot
- fn_out: text string with full path of output file name
- auto_vlims: a boolean governing how color limits are set
OUTPUT: either a screen image or a graphics file

"""

# imports

# The calling function, pan_plot.py, has already put alpha on the
# path so it is on the path here too.
import Lfun
Ldir = Lfun.Lstart()
if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
else: # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
    
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pickle
    
import zfun
import zrfun

from importlib import reload
import pfun; reload(pfun)
import pinfo; reload(pinfo)

def P_basic(in_dict):

    # START
    fig = plt.figure(figsize=(16,12)) # or pinfo.figsize for default
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
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_basic_salish(in_dict):

    # START
    fig = plt.figure(figsize=(20,10)) # or pinfo.figsize for default
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['salt', 'temp']
    aa = [-124, -122, 47, 49]
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
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #def add_windstress_flower(ax, ds, t_scl=0.2, t_leglen=0.1, center=(.85,.25)):
                # ADD MEAN WINDSTRESS VECTOR
                # t_scl: scale windstress vector (smaller to get longer arrows)
                # t_leglen: # Pa for wind stress vector legend
            pfun.add_windstress_flower(ax, ds, t_scl=.6, t_leglen=0.1, center=(.25, .3))
        elif ii == 2:
            #add_velocity_vectors(ax, ds, fn, v_scl=3, v_leglen=0.5, nngrid=80, zlev='top', center=(.7,.05)):
                # v_scl: scale velocity vector (smaller to get longer arrows)
                # v_leglen: m/s for velocity vector legend
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], v_scl=10, v_leglen=1.5, center=(.1, .1))
        ii += 1
    fig.tight_layout()
    
    # FINISH
    ds.close()
    if len(in_dict['fn_out']) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_phyt(in_dict):
    # a custom movie for LuAnne
    # START
    fig = plt.figure(figsize=pinfo.figsize)
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn_list = ['phytoplankton', 'phytoplankton']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn])
        pfun.add_coast(ax)
        aa = [-123.5, -122.1, 47, 49]
        pad = 0.005
        aap = [-123.5-pad, -122.1+pad, 47-pad, 49+pad]
        ax.plot([aa[0], aa[1], aa[1], aa[0], aa[0]], [aa[2], aa[2], aa[3], aa[3], aa[2]],
            '-m', linewidth=3)
        if ii == 1:
            fig.colorbar(cs)
            pfun.add_bathy_contours(ax, ds, txt=True)
            ax.axis(pfun.get_aa(ds))
        elif ii == 2:
            ax.axis(aap)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        elif ii == 2:
            ax.set_title('Puget Sound')
        ii += 1
        
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

def P_sect(in_dict):
    # plots a map and a section (distance, z)
    
    # START
    fig = plt.figure(figsize=(24,8))
    ds = nc.Dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'NO3'
    #vn = 'phytoplankton'
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
