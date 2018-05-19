"""
Functions for particle tracking.

Performance notes: in 3D it takes about:
 - cascadia1:  20 sec/day
 - cas3:      215 sec/day

"""
# setup
import numpy as np
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:sys.path.append(alp)
import zfun
import zrfun
import netCDF4 as nc4


def get_tracks(fn_list, plon0, plat0, pcs0, TR, trim_loc=False):
    
    # unpack items needed from TR
    
    if TR['rev']:
        dir_tag = 'reverse'
    else:
        dir_tag = 'forward'
    
    surface = not TR['3d']
    
    turb = TR['turb']
    ndiv = TR['ndiv']
    
    windage = TR['windage']
    
    # get basic info
    G = zrfun.get_basic_info(fn_list[0], only_G=True)
    maskr = np.ones_like(G['mask_rho']) # G['mask_rho'] = True in water
    maskr[G['mask_rho']==False] = 0 # maskr = 1 in water
    S = zrfun.get_basic_info(fn_list[0], only_S=True)

    # get time vector of history files
    NT = len(fn_list)
    rot = np.nan * np.ones(NT)
    counter = 0
    for fn in fn_list:
        ds = nc4.Dataset(fn)
        rot[counter] = ds.variables['ocean_time'][:].squeeze()
        counter += 1
        ds.close
    delta_t = rot[1] - rot[0] # second between saves

    # this is how we track backwards in time
    if dir_tag == 'reverse':
        delta_t = -delta_t
        fn_list = fn_list[::-1]

    # make vectors to feed to interpolant maker
    R = dict()
    R['rlonr'] = G['lon_rho'][0,:].squeeze()
    R['rlatr'] = G['lat_rho'][:,0].squeeze()
    R['rlonu'] = G['lon_u'][0,:].squeeze()
    R['rlatu'] = G['lat_u'][:,0].squeeze()
    R['rlonv'] = G['lon_v'][0,:].squeeze()
    R['rlatv'] = G['lat_v'][:,0].squeeze()
    R['rcsr'] = S['Cs_r'][:]
    R['rcsw'] = S['Cs_w'][:]
    
    # minimum grid sizes used when particles approach land boundaries
    dxg = np.diff(R['rlonr']).min()
    dyg = np.diff(R['rlatr']).min()
    
    # critical value of the landmask to decide if we are too close to land
    # closer to 1 excludes more particles
    # closer to 0 excludes more but they get zero velocities often
    maskr_crit = 0.8 # (maskr = 1 in water, 0 on land) [0.8 seems good]
    
    # these lists are used internally to get other variables as needed
    # (they must all be present in the history files)
    vn_list_vel = ['u','v','w']
    vn_list_zh = ['zeta','h']
    vn_list_wind = ['Uwind','Vwind']
    vn_list_other = ['salt', 'temp'] + vn_list_zh + vn_list_vel
    if windage == True:
        vn_list_other = vn_list_other + vn_list_wind
    # plist main is what ends up written to output
    plist_main = ['lon', 'lat', 'cs', 'ot', 'z'] + vn_list_other
    
    # Step through times.
    #
    counter = 0
    # rot is a list of all the ocean times (sec) in the current file list(e.g. 25)
    # and pot is a single one of these
    for pot in rot[:-1]:

        # get time indices of surrounding (in time) history files
        it0, it1, frt = zfun.get_interpolant(
                np.array(pot), rot, extrap_nan=False)
        # get Datasets
        ds0 = nc4.Dataset(fn_list[it0[0]])
        ds1 = nc4.Dataset(fn_list[it1[0]])

        if counter == 0:
            if trim_loc == True:
                # remove points on land
                pmask = zfun.interp_scattered_on_plaid(plon0, plat0,
                    R['rlonr'], R['rlatr'], maskr, exnan=True)
                # keep only points with pmask >= maskr_crit
                pcond = pmask >= maskr_crit
                plon = plon0[pcond]
                plat = plat0[pcond]
                pcs = pcs0[pcond]
            else:
                plon = plon0.copy()
                plat = plat0.copy()
                pcs = pcs0.copy()
            # create result arrays
            NP = len(plon)
            P = dict()
            for vn in plist_main:
                P[vn] = np.nan * np.ones((NT,NP))

            # write positions to the results arrays
            P['lon'][it0,:] = plon
            P['lat'][it0,:] = plat
            if surface == True:
                pcs[:] = S['Cs_r'][-1]
            P['cs'][it0,:] = pcs
            P = get_properties(vn_list_other, ds0, it0, P, plon, plat, pcs, R, surface)
            
        delt = delta_t/ndiv
        for nd in range(ndiv):
            
            # save the IC to use for out-of-bounds backup
            plonC = plon.copy()
            platC = plat.copy()
            #pcsC = pcs.copy()

            fr0 = nd/ndiv
            fr1 = (nd + 1)/ndiv
            frmid = (fr0 + fr1)/2

            # RK4 integration
            V0, ZH0 = get_vel(vn_list_vel, vn_list_zh,
                                       ds0, ds1, plon, plat, pcs, R, fr0, surface)
            plon1, plat1, pcs1 = update_position(V0, ZH0, S, delt/2,
                                                 plon, plat, pcs, surface)
            V1, ZH1 = get_vel(vn_list_vel, vn_list_zh,
                                       ds0, ds1, plon1, plat1, pcs1, R, frmid, surface)
            plon2, plat2, pcs2 = update_position(V1, ZH1, S, delt/2,
                                                 plon, plat, pcs, surface)
            V2, ZH2 = get_vel(vn_list_vel, vn_list_zh,
                                       ds0, ds1, plon2, plat2, pcs2, R, frmid, surface)
            plon3, plat3, pcs3 = update_position(V2, ZH2, S, delt,
                                                 plon, plat, pcs, surface)
            V3, ZH3 = get_vel(vn_list_vel, vn_list_zh,
                                       ds0, ds1, plon3, plat3, pcs3, R, fr1, surface)
            # add windage, calculated from the middle time
            if (surface == True) and (windage > 0):
                Vwind = get_wind(vn_list_wind, ds0, ds1, plon, plat, pcs, R, frmid, surface)
                Vwind3 = np.concatenate((windage*Vwind,np.zeros((NP,1))),axis=1)
            else:
                Vwind3 = np.zeros((NP,3))
            plon, plat, pcs = update_position((V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3,
                                              (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                              S, delt,
                                              plon, plat, pcs, surface)
                                              
            # If particles have moved onto land, move them back.
            # - first find those that ended up on land
            pmask = zfun.interp_scattered_on_plaid(plon, plat, R['rlonr'], R['rlatr'],
                maskr, exnan=True)
            pcond = pmask < maskr_crit # a Boolean mask
            # move those on land back to their staring point, with a random
            # perturbation of a grid cell to keep them from getting trapped
            if len(pcond) > 0:
                rix = np.random.randint(-1,1,len(plon))
                riy = np.random.randint(-1,1,len(plon))
                plon[pcond] = plonC[pcond] + rix[pcond]*dxg
                plat[pcond] = platC[pcond] + riy[pcond]*dyg
                # and check again for any stragglers that are still on land
                pmask = zfun.interp_scattered_on_plaid(plon, plat, R['rlonr'], R['rlatr'],
                    maskr, exnan=True)
                pcond = pmask < maskr_crit
                plon[pcond] = plonC[pcond]
                plat[pcond] = platC[pcond]
            
            # add turbulence in two distinct timesteps
            if turb == True:
                # pull values of VdAKs and add up to 3-dimensions
                VdAKs = get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frmid, surface)
                VdAKs3 = np.concatenate((np.zeros((NP,2)), VdAKs[:,np.newaxis]), axis=1)
                # update position with 1/2 of AKs gradient
                plon, plat, pcs = update_position(VdAKs3/2, (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                              S, delt, plon, plat, pcs, surface)
                # update position with rest of turbulence
                Vturb = get_turb(ds0, ds1, VdAKs, delta_t, plon, plat, pcs, R, frmid, surface)
                Vturb3 = np.concatenate((np.zeros((NP,2)), Vturb[:,np.newaxis]), axis=1)
                plon, plat, pcs = update_position(Vturb3, (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                              S, delt, plon, plat, pcs, surface)

        # write positions to the results arrays
        P['lon'][it1,:] = plon
        P['lat'][it1,:] = plat
        if surface == True:
            pcs[:] = S['Cs_r'][-1]
        P['cs'][it1,:] = pcs
        P = get_properties(vn_list_other, ds1, it1, P, plon, plat, pcs, R, surface)

        ds0.close()
        ds1.close()

        counter += 1

    # by doing this the points are going forward in time
    if dir_tag == 'reverse':
        for vn in plist_main:
            P[vn] = P[vn][::-1,:]

    # and save the time vector (seconds in whatever the model reports)
    P['ot'] = rot

    return P

def update_position(V, ZH, S, delta_t, plon, plat, pcs, surface):
    # find the new position
    Plon = plon.copy()
    Plat = plat.copy()
    Pcs = pcs.copy()
    dX_m = V*delta_t
    per_m = zfun.earth_rad(Plat)
    clat = np.cos(np.pi*Plat/180.)
    pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
    pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
    H = ZH.sum(axis=1)
    # NOTE 2018.05.18 I believe that the first column of ZH is the
    # zeta of all particles, and the second is h, so the H we calculate
    # above is the total water column thickness.  The particle z-positions
    # should then be ...
    pdz_s = dX_m[:,2]/H
    Plon += pdx_deg
    Plat += pdy_deg
    if surface == False:
        Pcs_orig = Pcs.copy()
        Pcs += pdz_s
        # enforce limits on cs
        mask = np.isnan(Pcs) # does this happen?
        if len(mask) > 0:
            Pcs[mask] = Pcs_orig[mask]
        pad = 1e-5 # stay a bit away from extreme values
        Pcs[Pcs < S['Cs_r'][0]+pad] = S['Cs_r'][0]+pad
        Pcs[Pcs > S['Cs_r'][-1]-pad] = S['Cs_r'][-1]-pad
    else:
        Pcs[:] = S['Cs_r'][-1]

    return Plon, Plat, Pcs

def get_vel(vn_list_vel, vn_list_zh, ds0, ds1, plon, plat, pcs, R, frac, surface):
    # get the velocity, zeta, and h at all points, at an arbitrary
    # time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1
    # 0 <= frac <= 1
    V0 = get_V(vn_list_vel, ds0, plon, plat, pcs, R, surface)
    V1 = get_V(vn_list_vel, ds1, plon, plat, pcs, R, surface)
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R, surface)
    ZH1 = get_V(vn_list_zh, ds1, plon, plat, pcs, R, surface)
    V = (1 - frac)*V0 + frac*V1
    ZH = (1 - frac)*ZH0 + frac*ZH1
    
    return V, ZH

def get_wind(vn_list_wind, ds0, ds1, plon, plat, pcs, R, frac, surface):
    # get the wind velocity at an arbitrary
    # time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1
    # 0 <= frac <= 1
    V0 = get_V(vn_list_wind, ds0, plon, plat, pcs, R, surface)
    V1 = get_V(vn_list_wind, ds1, plon, plat, pcs, R, surface)
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    V = (1 - frac)*V0 + frac*V1

    return V

def get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frac, surface):
    # create diffusivity gradient for turbulence calculation
    
    # first time step
    ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R, surface)
    ZH0[np.isnan(ZH0)] = 0
    dpcs0 = 1/(ZH0[:,0] + ZH0[:,1]) # change in pcs for a total of a 2m difference
    
    #     upper variables
    pcs0u = pcs-dpcs0
    pcs0u[pcs0u > S['Cs_r'][-1]] = S['Cs_r'][-1]
    AKs0u = get_V(['AKs',], ds0, plon, plat, pcs0u, R, surface) # diffusivity
    z0u = (pcs0u)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
    
    #     lower variables
    pcs0b = pcs+dpcs0
    pcs0b[pcs0b < S['Cs_r'][0]] = S['Cs_r'][0]
    AKs0b = get_V(['AKs',], ds0, plon, plat, pcs0b, R, surface) # diffusivity
    z0b = (pcs0b)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
    
    #     combine at midpoint
    V0 = (AKs0u-AKs0b).squeeze()/(z0u-z0b)
    
    # second time step
    ZH1 = get_V(vn_list_zh, ds1, plon, plat, pcs, R, surface)
    ZH1[np.isnan(ZH1)] = 0
    dpcs1 = 1/(ZH1[:,0] + ZH1[:,1]) # change in pcs for 1m difference
    
    #     upper variables
    pcs1u = pcs-dpcs1
    pcs1u[pcs0u > S['Cs_r'][-1]] = S['Cs_r'][-1]
    AKs1u = get_V(['AKs',], ds1, plon, plat, pcs1u, R, surface) # diffusivity
    z1u = (pcs1u)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
    
    #     lower variables
    pcs1b = pcs+dpcs0
    pcs1b[pcs1b < S['Cs_r'][0]] = S['Cs_r'][0]
    AKs1b = get_V(['AKs',], ds1, plon, plat, pcs0b, R, surface) # diffusivity
    z1b = (pcs1b)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
    
    #     combine at midpoint
    V1 = (AKs1u-AKs1b).squeeze()/(z1u-z1b)
    
    # average of timesteps
    V = (1-frac)*V0 + frac*V1
    
    return V
    
def get_turb(ds0, ds1, dAKs, delta_t, plon, plat, pcs, R, frac, surface):
    # get the vertical turbulence correction components
    
    # getting diffusivity
    V0 = get_V(['AKs',], ds0, plon, plat, pcs, R, surface).squeeze()
    V1 = get_V(['AKs',], ds1, plon, plat, pcs, R, surface).squeeze()
    # replace nans
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    # create weighted average diffusivity
    Vave = (1 - frac)*V0 + frac*V1
    
    # turbulence calculation from Banas, MacCready, and Hickey (2009)
    # w_turbulence = rand*sqrt(2k/delta(t)) + delta(k)/delta(z)
    # rand = random array with normal distribution
    rand = np.random.standard_normal(len(V0))
    # only using half of gradient, first half already added
    V = rand*np.sqrt(2*Vave/delta_t) + dAKs/2
    
    return V

def get_properties(vn_list_other, ds, it, P, plon, plat, pcs, R, surface):
    # find properties at a position
    OTH = get_V(vn_list_other, ds, plon, plat, pcs, R, surface)
    for vn in vn_list_other:
        P[vn][it,:] = OTH[:,vn_list_other.index(vn)]
    this_zeta = OTH[:, vn_list_other.index('zeta')]
    this_h = OTH[:, vn_list_other.index('h')]
    full_depth = this_zeta + this_h
    P['z'][it,:] = pcs * full_depth # is this correct??

    return P

def get_V(vn_list, ds, plon, plat, pcs, R, surface):

    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages

    # get interpolant arrays
    i0lon_d = dict()
    i1lon_d = dict()
    frlon_d = dict()
    i0lat_d = dict()
    i1lat_d = dict()
    frlat_d = dict()
    
    for gg in ['r', 'u', 'v']:
        exn = True # nan-out particles that leave domain
        i0lon_d[gg], i1lon_d[gg], frlon_d[gg] = zfun.get_interpolant(
                plon, R['rlon'+gg], extrap_nan=exn)
        i0lat_d[gg], i1lat_d[gg], frlat_d[gg] = zfun.get_interpolant(
                plat, R['rlat'+gg], extrap_nan=exn)
    i0csr, i1csr, frcsr = zfun.get_interpolant(pcs, R['rcsr'], extrap_nan=exn)
    i0csw, i1csw, frcsw = zfun.get_interpolant(pcs, R['rcsw'], extrap_nan=exn)
    
    NV = len(vn_list)
    NP = len(plon)
    
    # get interpolated values
    V = np.nan * np.ones((NP,NV))
    vcount = 0
    for vn in vn_list:
        if vn in ['w']:
            i0cs = i0csw
            i1cs = i1csw
            frcs = frcsw
        else:
            i0cs = i0csr
            i1cs = i1csr
            frcs = frcsr
        if vn == 'u':
            gg = 'u'
        elif vn == 'v':
            gg = 'v'
        else:
            gg = 'r'
        i0lat = i0lat_d[gg]
        i1lat = i1lat_d[gg]
        frlat = frlat_d[gg]
        i0lon = i0lon_d[gg]
        i1lon = i1lon_d[gg]
        frlon = frlon_d[gg]
        # get the data field and put nan's in masked points
        # (later we do more massaging)
        if vn in ['salt','temp','u','v','w'] and surface==True:
            v0 = ds[vn][0, -1, :, :].squeeze()
        else:
            v0 = ds[vn][:].squeeze()
        try:
            vv = v0.data
            vv[v0.mask] = np.nan
        except AttributeError:
            # it is not a masked array
            vv = v0
            
        if vn in ['salt','temp','AKs','u','v','w'] and surface==False:
            # Get just the values around our particle positions.
            # each row in VV corresponds to a "box" around a point
            #
            VV = np.nan* np.ones((NP, 8))
            # using "fancy indexing" in which the three indexing arrays
            # are the same length and the resulting vector has the same length
            # (works for numpy arrays, not for NetCDF)
            VV[:,0] = vv[i0cs, i0lat, i0lon]
            VV[:,1] = vv[i0cs, i0lat, i1lon]
            VV[:,2] = vv[i0cs, i1lat, i0lon]
            VV[:,3] = vv[i0cs, i1lat, i1lon]
            VV[:,4] = vv[i1cs, i0lat, i0lon]
            VV[:,5] = vv[i1cs, i0lat, i1lon]
            VV[:,6] = vv[i1cs, i1lat, i0lon]
            VV[:,7] = vv[i1cs, i1lat, i1lon]
            # Work on edge values.  If all in a box are masked
            # then that row will be 8 nan's, and also:
            if vn in ['u', 'v', 'w', 'AKs']:
                # set all velocities and AKs to zero if any in the box are nan
                mask = np.isnan(VV)
                VV[mask] = 0
            elif vn in ['salt','temp']:
                # set all tracers to their average if any in the box are masked
                newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,8))
                mask = np.isnan(VV)
                VV[mask] = newval[mask]
            # now do the interpolation in each box
            vl = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
            vu = ( (1-frlat)*((1-frlon)*VV[:,4] + frlon*VV[:,5])
                + frlat*((1-frlon)*VV[:,6] + frlon*VV[:,7]) )
            v = (1-frcs)*vl + frcs*vu
        elif vn in ['salt','temp','u','v','w'] and surface==True:
            VV = np.nan* np.ones((NP, 4))
            VV[:,0] = vv[i0lat, i0lon]
            VV[:,1] = vv[i0lat, i1lon]
            VV[:,2] = vv[i1lat, i0lon]
            VV[:,3] = vv[i1lat, i1lon]
            # Work on edge values.  If all in a box are masked
            # then that row will be nan's, and also:
            if vn in ['u', 'v', 'w']:
                # set all velocities to zero if any in the box are masked
                mask = np.isnan(VV)
                VV[mask] = 0
            elif vn in ['salt','temp']:
                # set all tracers to their average if any in the box are masked
                newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
                mask = np.isnan(VV)
                VV[mask] = newval[mask]
            v = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
        elif vn in ['zeta','Uwind','Vwind', 'h']:
            VV = np.nan* np.ones((NP, 4))
            VV[:,0] = vv[i0lat, i0lon]
            VV[:,1] = vv[i0lat, i1lon]
            VV[:,2] = vv[i1lat, i0lon]
            VV[:,3] = vv[i1lat, i1lon]
            newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
            mask = np.isnan(VV)
            VV[mask] = newval[mask]
            v = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
        V[:,vcount] = v
        vcount += 1
            
    return V

def get_fn_list(idt, Ldir):
    # LiveOcean version, for 1 day only.
    # Assumes we have history files 1-25, corresponding to hours 0-24.
    fn_list = []
    dd = idt.strftime('%Y.%m.%d')
    indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
            '/f' + dd + '/')
    for hh in range(1,26):
        hhhh = ('0000' + str(hh))[-4:]
        fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
    return fn_list