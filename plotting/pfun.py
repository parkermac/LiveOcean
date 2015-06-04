"""
Module of plotting functions.
"""
def roms_basic(fn, alp, fn_coast='', show_plot=True, save_plot=False,
    fn_out='test.png'):    
    # This creates, and optionally saves, a basic plot of surface fields
    # from a ROMS history file.
       
    # INPUT   
    # fn: text string with the full path name of the history file to plot   
    # fn_coast: text string with the full path to a coastline data file    
    # alp: text string of the path to where the modules zfun and matfun are    
    # show_plot: Boolean, set to True to show plot on screen    
    # save_plot: Boolean, set to True to save plot to a file    
    # fn_out: text string with full path of output file name
    
    # OUTPUT
    # either a screen image or a graphics file
    
    # scaling choices
    #sname = 'south_sound'
    sname = 'cascadia'
    if sname == 'south_sound':
        salt_lims = (27,31)
        temp_lims = (10,22)  
        v_scl = 30 # scale velocity vector (smaller to get longer arrows)
        v_leglen = 1 # m/s for velocity vector legend
        t_scl = .2 # scale windstress vector (smaller to get longer arrows)
        t_leglen = 0.1 # Pa for wind stress vector legend
    elif sname == 'cascadia':
        salt_lims = (28, 34)
        temp_lims = (6, 18)    
        v_scl = 3 # scale velocity vector (smaller to get longer arrows)
        v_leglen = 0.5 # m/s for velocity vector legend
        t_scl = .2 # scale windstress vector (smaller to get longer arrows)
        t_leglen = 0.1 # Pa for wind stress vector legend
        
    # setup
    import sys
    if alp not in sys.path:
        sys.path.append(alp)    
    import zfun; reload(zfun) # utility functions
    import matfun; reload(matfun) # functions for working with mat files
    import numpy as np

    # GET DATA
    G, S, T = zfun.get_basic_info(fn)
    import netCDF4 as nc   
    ds = nc.Dataset(fn,'r')
    h = G['h']
    salt = ds.variables['salt'][0, -1, :, :].squeeze()
    temp = ds.variables['temp'][0, -1, :, :].squeeze()
    u = ds.variables['u'][0, -1, :, :].squeeze()
    v = ds.variables['v'][0, -1, :, :].squeeze()  
    taux = ds.variables['sustr'][:].squeeze()
    tauy = ds.variables['svstr'][:].squeeze()  
    ds.close()     
    lonp = G['lon_psi']
    latp = G['lat_psi']
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]   
    depth_levs = [100, 200, 500, 1000, 2000, 3000]
       
    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)
        
    # PLOTTING       
    import matplotlib.pyplot as plt
    plt.close()
    #fig = plt.figure(figsize=(14, 8))
    fig = plt.figure(figsize=(14, 8))
    
    # 1. surface salinity    
    ax = fig.add_subplot(121)
    cmap = plt.get_cmap(name='jet')    
    # pcolormesh is a fast way to make pcolor plots.  By handing it lon_psi
    # and lat_psi and any rho_grid_layer[1:-1,1:-1] the coordinate arrays
    # are one bigger in size (in both dimensions) than the data array, meaning
    # they define the data corners.  The result is that the colored tiles are
    # centered exactly where they should be (at their rho-grid location).
    cs = ax.pcolormesh(lonp, latp, salt[1:-1,1:-1],
        vmin=salt_lims[0], vmax=salt_lims[1],  cmap = cmap)    
    # bathymetry contours
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')        
    # add coastline
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)        
    # extras
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Surface Salinity')
    fig.colorbar(cs)
    # put info on plot
    ax.text(.95, .07, T['tm'].strftime('%Y-%m-%d'),
        horizontalalignment='right', transform=ax.transAxes)
    ax.text(.95, .04, T['tm'].strftime('%H:%M') + ' UTC',
        horizontalalignment='right', transform=ax.transAxes)
        
    # ADD MEAN WINDSTRESS VECTOR
    tauxm = taux.mean()
    tauym = tauy.mean()
    ax.quiver([.85, .85] , [.25, .25], [tauxm, tauxm], [tauym, tauym],
        units='y', scale=t_scl, scale_units='y', color='k', transform=ax.transAxes)
    tt = 1./np.sqrt(2)
    t_alpha = 0.3
    ax.quiver([.85, .85] , [.25, .25],
        t_leglen*np.array([0,tt,1,tt,0,-tt,-1,-tt]),
        t_leglen*np.array([1,tt,0,-tt,-1,-tt,0,tt]),
        units='y', scale=t_scl, scale_units='y', color='k', alpha=t_alpha, transform=ax.transAxes)
    ax.text(.85, .12,'Windstress',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes)
    ax.text(.85, .15, str(t_leglen) + ' Pa',
        horizontalalignment='center', alpha=t_alpha, transform=ax.transAxes)
    
    # 2. surface temperature    
    ax = fig.add_subplot(122)
    cmap = plt.get_cmap(name='rainbow')
    cs = ax.pcolormesh(lonp, latp, temp[1:-1,1:-1],
        vmin=temp_lims[0], vmax=temp_lims[1],  cmap = cmap)        
    # bathymetry
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')        
    # add coastline
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)        
    # extras
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Surface Temperature $^{\circ}C$')
    fig.colorbar(cs)
        
    # ADD VELOCITY VECTORS    
    # set masked values to 0
    ud = u.data; ud[G['mask_u']==False] = 0
    vd = v.data; vd[G['mask_v']==False] = 0    
    # create interpolant
    import scipy.interpolate as intp
    ui = intp.interp2d(G['lon_u'][0, :], G['lat_u'][:, 0], ud)
    vi = intp.interp2d(G['lon_v'][0, :], G['lat_v'][:, 0], vd)    
    # create regular grid
    import numpy as np
    aaa = ax.axis()
    daax = aaa[1] - aaa[0]
    daay = aaa[3] - aaa[2]
    axrat = np.cos(np.deg2rad(aaa[2])) * daax / daay
    nngrid = 80
    x = np.linspace(aaa[0], aaa[1], round(nngrid * axrat))
    y = np.linspace(aaa[2], aaa[3], nngrid)
    xx, yy = np.meshgrid(x, y)    
    # interpolate to regular grid
    uu = ui(x, y)
    vv = vi(x, y) 
    mask = uu != 0    
    # plot velocity vectors  
    ax.quiver(xx[mask], yy[mask], uu[mask], vv[mask],
        units='y', scale=v_scl, scale_units='y', color='k')
    ax.quiver([.7, .7] , [.05, .05], [v_leglen, v_leglen], [v_leglen, v_leglen],
        units='y', scale=v_scl, scale_units='y', color='k', transform=ax.transAxes)
    ax.text(.75, .05, str(v_leglen) + ' $ms^{-1}$',
        horizontalalignment='left', transform=ax.transAxes)
    # END OF VELOCITY VECTOR SECTION
         
    if show_plot==True:
        plt.show()   
    if save_plot==True:
        plt.savefig(fn_out)

def roms_layer(fn, alp, fn_coast='', show_plot=True, save_plot=False,
    fn_out='test.png'):
    # plots fields on a specified depth level
           
    # setup
    import sys
    if alp not in sys.path:
        sys.path.append(alp)    
    import zfun; reload(zfun) # utility functions
    import matfun; reload(matfun) # functions for working with mat files
    
    zlev = -500. # z (m) of layer to plot
    which_var = 'salt'
    
    # IMPORTS
    import matplotlib.pyplot as plt
    import numpy as np
    
    # GET DATA    
    G, S, T = zfun.get_basic_info(fn)
    h = G['h']
    M = G['M']
    L = G['L']
    lonp = G['lon_psi']
    latp = G['lat_psi']
    lonr = G['lon_rho']
    latr = G['lat_rho']
    mask = G['mask_rho']   
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    depth_levs = [100, 200, 500, 1000, 2000, 3000]
       
    # make full z field (append top and bottom)
    import netCDF4 as nc   
    ds = nc.Dataset(fn,'r')
    # I set zeta to zero below because that is better when making layers
    # (e.g. "zlev = -1." will be more like 1 m below the surface than if
    # we had used the actual zeta)
    zeta = 0 * ds.variables['zeta'][:].squeeze()
    zr_mid = zfun.get_z(h, zeta, S, only_rho=True)
    zr_bot = -h.reshape(1, M, L).copy()
    zr_top = zeta.reshape(1, M, L).copy()
    zr = zfun.make_full((zr_bot, zr_mid, zr_top))
    
    # get the data field and its units - if any
    fld_mid = ds.variables[which_var][:].squeeze()
    try:
        fld_units = str(ds.variables[which_var].units)
    except:
        fld_units = ''
    # make full data field (extrapolate top and bottom)
    fld = zfun.make_full((fld_mid,))
    
    ds.close()
    
    # convert units if needed
    if which_var == 'oxygen':
        fld = 1000.0 * fld / 44660
        fld_units = 'mL/L'
    
    # INTERPOLATION
    which_z = zlev * np.ones(1)
    lay = zfun.get_layer(fld, zr, which_z)
    
    # get the surface field to plot
    print str(fld.shape)
    lay_top = fld[-1].copy() # fld[-1] is shorthand for fld[-1,:,:]
    
    # mask out the bio fields if we are inside the Salish Sea, because
    # we know that those equations were not activated there
    if which_var in ['NO3', 'phytoplankton', 'zooplankton', 'detritus',
    'Ldetritus', 'oxygen']:
        # path made using sal = plt.ginput(0), RETURN to end,
        # and then adjusted by hand
        sal = [(-122.0134210053738, 50.1),
        (-122.0134210053738, 46.8),
        (-123.33921283879502, 46.9),
        (-124.57460977448297, 48.249154246315108),
        (-125.26763732377132, 50.1)]
        from matplotlib.path import Path
        p = Path(sal)
        pin = p.contains_points(np.array(zip(lonr.flatten(), latr.flatten())))
        pin = pin.reshape(lonr.shape) # a Boolean array
        #
        lay_top[pin == True] = np.nan
        lay[pin == True] = np.nan
    # make masking on layers to be plotted (to work with pcolormesh)
    import numpy.ma as ma
    lay_top[mask == False] = np.nan
    laym_top = ma.masked_where(np.isnan(lay_top), lay_top)
    # also mask the horizontal slice while we are at it
    lay[mask == False] = np.nan
    laym = ma.masked_where(np.isnan(lay), lay)
        
    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)
               
    # PLOTTING    
    plt.close()
    fig = plt.figure(figsize=(16, 10))
    cmap = plt.get_cmap(name='rainbow')
    
    # color limits
    if which_var == 'salt':
        clim_top = (29,34)
        clim = (34,35)
    else:
        clim_top = ( np.floor(laym_top.min()), np.ceil(laym_top.max()) )
        clim= ( np.floor(laym.min()), np.ceil(laym.max()) )
     
    # (1) map of surface field
    ax = fig.add_subplot(121)    
    cs = ax.pcolormesh(lonp, latp, laym_top[1:-1,1:-1], cmap=cmap) 
    cs.set_clim(clim_top)    
    # bathymetry
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')        
    # add coastline
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)    
    # limits, scaling, and labels
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('(a) Surface ' + which_var + ' ' + fld_units)
    fig.colorbar(cs)    
    # put info on plot
    ax.text(.95, .1, T['tm'].strftime('%Y-%m-%d'),
        horizontalalignment='right', transform=ax.transAxes)
    ax.text(.95, .05, T['tm'].strftime('%H:%M'),
        horizontalalignment='right', transform=ax.transAxes)
       
    # (2) map of the horizontal layer
    ax = fig.add_subplot(122)    
    # plot the field using pcolormesh
    cs = ax.pcolormesh(lonp, latp, laym[1:-1,1:-1],  cmap=cmap) 
    cs.set_clim(clim)   
    # bathymetry
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')        
    # add coastline
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)    
    # limits, scaling, and labels
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_title('(b) ' + which_var + ' at ' + str(zlev) + ' m')
    fig.colorbar(cs)
    
    if show_plot==True:
        plt.show()   
    if save_plot==True:
        plt.savefig(fn_out)
        
def roms_sect(fn, alp, fn_coast='', show_plot=True, save_plot=False,
    fn_out='test.png'):
    # plots a section (distance, z)
    
    # setup
    import os; import sys
    alp = os.path.abspath('../alpha')
    if alp not in sys.path:
        sys.path.append(alp)
    import Lfun; reload(Lfun)
    Ldir = Lfun.Lstart(alp)
    
    import matplotlib.pyplot as plt
    import numpy as np   
    import zfun; reload(zfun) # plotting functions
    import matfun; reload(matfun) # functions for working with mat files
          
    # GET DATA
    G, S, T = zfun.get_basic_info(fn)
    import netCDF4 as nc   
    ds = nc.Dataset(fn,'r')
    h = G['h']
    zeta = ds.variables['zeta'][:].squeeze()
    zr = zfun.get_z(h, zeta, S, only_rho=True)
    
    varname = 'temp'
    
    sectvar = ds.variables[varname][:].squeeze()
    clims = {'salt':(28,35),'temp':(4,14)}
    
    L = G['L']
    M = G['M']
    N = S['N']
    
    lon = G['lon_rho']
    lat = G['lat_rho']
    mask = G['mask_rho']
    maskr = mask.reshape(1, M, L).copy()
    mask3 = np.tile(maskr, [N, 1, 1])
    zbot = -h # don't need .copy() because of the minus operation
    
    # make sure fields are masked
    zeta[mask==False] = np.nan
    zbot[mask==False] = np.nan
    sectvar[mask3==False] = np.nan
    
    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)
    
    # CREATE THE SECTION    
    # which track to use for the section
    if False:
        tracks_path = Ldir['data'] + 'tracks/'
        which_track = 'lat47Track'
        #which_track = 'jdf2psTrack'
        # get the track to interpolate onto
        mat = matfun.loadmat(tracks_path + which_track + '.mat')
        x = mat['x']
        y = mat['y']
    else:   
        # create track by hand
        which_track = 'By hand'
        x = np.linspace(lon.min(), -124, 500)
        y = 48.5 * np.ones(x.shape)
    # create dist   
    earth_rad = 6371e3 # m
    xrad = np.pi * x /180
    yrad = np.pi * y / 180
    dx = earth_rad * np.cos(yrad[1:]) * np.diff(xrad)
    dy = earth_rad * np.diff(yrad)
    ddist = np.sqrt(dx**2 + dy**2)
    dist = np.zeros(len(x)) 
    dist[1:] = ddist.cumsum()/1000 # km        
    # find the index of zero
    it0 = zfun.get_interpolant(np.zeros(1), dist)
    idist0 = int(it0[0][0]) + int((it0[0][2].round()))    
    distr = dist.reshape(1, len(dist)).copy()
    dista = np.tile(distr, [N, 1]) # array    
    # pack fields to process in dicts
    d2 = dict()
    d2['zbot'] = zbot
    d2['zeta'] = zeta
    d2['lon'] = lon
    d2['lat'] = lat
    d3 = dict()
    d3['zr'] = zr
    d3['sectvar'] = sectvar    
    # get vectors describing the (plaid) grid
    xx = lon[1,:]
    yy = lat[:,1]
    xit = zfun.get_interpolant(x, xx)
    yit = zfun.get_interpolant(y, yy)    
    # and prepare them to do the bilinear interpolation
    xita = np.array(xit)
    yita = np.array(yit)
    col0 = xita[:, 0].astype(int)
    col1 = xita[:, 1].astype(int)
    colf = xita[:, 2]
    colff = 1 - colf
    row0 = yita[:, 0].astype(int)
    row1 = yita[:, 1].astype(int)
    rowf = yita[:, 2]
    rowff = 1 - rowf    
    # now actually do the interpolation    
    # 2-D fields
    v2 = dict()
    for fname in d2.keys():
        fld = d2[fname]
        fldi = (rowff*(colff*fld[row0, col0] + colf*fld[row0, col1]) 
        + rowf*(colff*fld[row1, col0] + colf*fld[row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldi = fldi.data # just the data, not the mask
        v2[fname] = fldi    
    # 3-D fields
    v3 = dict()
    for fname in d3.keys():
        fld = d3[fname]
        fldi = (rowff*(colff*fld[:, row0, col0] + colf*fld[:, row0, col1]) 
        + rowf*(colff*fld[:, row1, col0] + colf*fld[:, row1, col1]))
        if type(fldi) == np.ma.core.MaskedArray:
            fldi = fldi.data # just the data, not the mask
        v3[fname] = fldi        
    v3['dist'] = dista # distance in km    
    # make "full" fields by padding top and bottom
    nana = np.nan * np.ones((N + 2, len(dist))) # blank array
    v3['zrf'] = nana.copy()
    v3['zrf'][0,:] = v2['zbot']
    v3['zrf'][1:-1,:] = v3['zr']
    v3['zrf'][-1,:] = v2['zeta']
    #
    v3['sectvarf'] = nana.copy()
    v3['sectvarf'][0,:] = v3['sectvar'][0,:]
    v3['sectvarf'][1:-1,:] = v3['sectvar']
    v3['sectvarf'][-1,:] = v3['sectvar'][-1,:]
    #
    v3['distf'] = nana.copy()
    v3['distf'][0,:] = v3['dist'][0,:]
    v3['distf'][1:-1,:] = v3['dist']
    v3['distf'][-1,:] = v3['dist'][-1,:]    
    # NOTE: should make this a function
    # (but note that it is specific to each variable)
    
    # PLOTTING    
    plt.close()
    fig = plt.figure(figsize=(20,10))
    
    # map
    ax = fig.add_subplot(131)
    cmap = plt.get_cmap(name='Set3')
    ax.contourf(lon, lat, zbot, 50, vmin=-3000, vmax=0, cmap=cmap)
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    ax.contour(lon, lat, zbot, [-2000, -1000, -500, -200, -100], colors='g')
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=10, markerfacecolor='w',
    markeredgecolor='r', markeredgewidth=2)
    # test the interpolation routine
    #ax.plot(v2['lon'], v2['lat'], '*k') # RESULT: it works
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)    
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Bathymetry and Section Track')
    
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.set_ylim(v3['zr'].min(), 5)
    ax.set_ylim(-1000, 5)
    cmap = plt.get_cmap(name='rainbow')
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], v3['sectvarf'],  cmap=cmap) 
    cs.set_clim(clims[varname])   
    fig.colorbar(cs)
    cs = ax.contour(v3['distf'], v3['zrf'], v3['sectvarf'],
        np.arange(0, 40, .25), colors='k')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title(which_track + ' ' + varname)
    
    if show_plot==True:
        plt.show()   
    if save_plot==True:
        plt.savefig(fn_out)
        
def nest_plot(fn, alp, fn_coast='', show_plot=True, save_plot=False,
    fn_out='test.png'):
    # plots a field nested inside the corresponding HYCOM field
    
    # setup
    import os; import sys
    alp = os.path.abspath('../alpha')
    if alp not in sys.path:
        sys.path.append(alp)
    import Lfun; reload(Lfun)
    Ldir = Lfun.Lstart(alp)
    import zfun; reload(zfun)
    import matfun; reload(matfun)
    import netCDF4 as nc
    import numpy as np
    
    which_var = 'salt'
    # color limits
    if which_var == 'temp':
        v0 = 6; v1 = 12
    elif which_var == 'salt':
        v0 = 33.8; v1 = 34.2;

    # valid depth choices (start,stop,step): (0,12,2), (15,50,5), (60,100,10),
    # 125, (150,400,50), (500,1000,100), 1250, (1500,3000,500), 4000, 5000
    zlev = -300.
    
    depth_levs = [100, 200, 500, 1000, 2000, 3000] # bathy plotting
       
    # HYCOM
    hyvar_dict = {'temp':'t3d', 'salt':'s3d'}
    hvar = hyvar_dict[which_var]
    f_string = fn.split('/')[-2]
    fnhy = (Ldir['LOo'] + Ldir['gtag'] + '/' + f_string +
        '/ocn/Data/' + hvar + '.nc')
    ds = nc.Dataset(fnhy)
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    z = ds.variables['z'][:]
    nz = list(z).index(zlev)
    laym_hy = ds.variables[hvar + '_filt'][0,nz,:,:].squeeze()
    ds.close()
    
    # CLM
    #fnclm = (Ldir['LOo'] + Ldir['gtag'] + '/' + f_string +
    #    '/ocn/ocean_clm.nc') 
    #ds = nc.Dataset(fnclm)
    #zfun.ncd(ds)
    #ds.close()
    
    # HISTORY
    ds = nc.Dataset(fn)
    # get info
    G, S, T = zfun.get_basic_info(fn)
    h = G['h']
    M = G['M']
    L = G['L']
    lonp = G['lon_psi']
    latp = G['lat_psi']
    mask = G['mask_rho']
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    
    # get a horizontal layer
    zeta = 0 * h
    zr_mid = zfun.get_z(h, zeta, S, only_rho=True)
    zr_bot = -h.reshape(1, M, L).copy()
    zr_top = zeta.reshape(1, M, L).copy()
    zr = zfun.make_full((zr_bot, zr_mid, zr_top))
    # get the data field
    fld_mid = ds.variables[which_var][:].squeeze()
    # make full data field (extrapolate top and bottom)
    fld = zfun.make_full((fld_mid,))
    # INTERPOLATION
    #import numpy as np
    which_z = zlev * np.ones(1)
    lay = zfun.get_layer(fld, zr, which_z)
    # make masking on layers to be plotted (to work with pcolormesh)
    import numpy.ma as ma
    lay[mask == False] = np.nan
    laym = ma.masked_where(np.isnan(lay), lay)
    
    # get coastline
    if len(fn_coast) != 0:
        # get the coastline
        cmat = matfun.loadmat(fn_coast)
        
    # PLOTTING   
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(name='rainbow')
    plt.close()
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 10), squeeze=False)
        
    ax = axes[0,0]
    ax.set_title('HYCOM: z=' + str(zlev) + 'm, var=' + which_var)
    cs = ax.pcolormesh(lon, lat, laym_hy, vmin=v0, vmax=v1, cmap = cmap)
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    ax.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    aa = ax.axis()
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    fig.colorbar(cs, ax=ax)
    
    ax = axes[0,1]
    ax.set_title('HYCOM + ROMS')
    cs = ax.pcolormesh(lon, lat, laym_hy, vmin=v0, vmax=v1, cmap = cmap)
    ax.pcolormesh(lonp, latp, laym[1:-1,1:-1], vmin=v0, vmax=v1, cmap = cmap)
    ax.contour(G['lon_rho'], G['lat_rho'], h, depth_levs, colors='g')
    if len(fn_coast) != 0:
        ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    ax.axis(aa)
    zfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    fig.colorbar(cs, ax=ax)

    if show_plot==True:
        plt.show()   
    if save_plot==True:
        plt.savefig(fn_out)

