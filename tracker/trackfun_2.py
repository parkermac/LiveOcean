"""
Functions for particle tracking.

I am in the process of moving this code over to using nearest-neighbor
interpolation vor everything.

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
from scipy.spatial import cKDTree
import pickle

from time import time

# Shared Constants
#
# criterion for deciding if particles are on land
maskr_crit = 0.8 # (maskr = 1 in water, 0 on land) [0.8 seems good]

# NEW CODE for nearest neighbor interpolation
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'
fn = Ldir['roms'] + '/output/' + Ldir['gtagex'] + '/f2017.07.04/ocean_his_0001.nc'
G, S, T = zrfun.get_basic_info(fn)
Maskr = G['mask_rho'] # True over water
Masku = G['mask_u'] # True over water
Maskv = G['mask_v'] # True over water
Maskr3 = np.tile(G['mask_rho'].reshape(1,G['M'],G['L']),[S['N'],1,1])
Masku3 = np.tile(G['mask_u'].reshape(1,G['M'],G['L']-1),[S['N'],1,1])
Maskv3 = np.tile(G['mask_v'].reshape(1,G['M']-1,G['L']),[S['N'],1,1])
Maskw3 = np.tile(G['mask_rho'].reshape(1,G['M'],G['L']),[S['N']+1,1,1])
# load trees
tree_dir = Ldir['LOo'] + 'tracker_trees/cas6/'
xyT_rho = pickle.load(open(tree_dir + 'xyT_rho.p', 'rb'))
xyT_rho_un = pickle.load(open(tree_dir + 'xyT_rho_un.p', 'rb'))
xyT_u = pickle.load(open(tree_dir + 'xyT_u.p', 'rb'))
xyT_v = pickle.load(open(tree_dir + 'xyT_v.p', 'rb'))
xyzT_rho = pickle.load(open(tree_dir + 'xyzT_rho.p', 'rb'))
xyzT_u = pickle.load(open(tree_dir + 'xyzT_u.p', 'rb'))
xyzT_v = pickle.load(open(tree_dir + 'xyzT_v.p', 'rb'))
xyzT_w = pickle.load(open(tree_dir + 'xyzT_w.p', 'rb'))

def get_tracks(fn_list, plon0, plat0, pcs0, TR, trim_loc=False):
    """
    This is the main function doing the particle tracking.
    """
    
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
    maskr = maskr.flatten()
    S = zrfun.get_basic_info(fn_list[0], only_S=True)
    
    # minimum grid sizes
    dxg = np.diff(G['lon_rho'][0,:]).min()
    dyg = np.diff(G['lat_rho'][:,0]).min()

    # get time vector of history files
    NT = len(fn_list)
    rot = np.nan * np.ones(NT)
    counter = 0
    for fn in fn_list:
        ds = nc4.Dataset(fn)
        rot[counter] = ds.variables['ocean_time'][:].squeeze()
        counter += 1
        ds.close
    delta_t_his = rot[1] - rot[0] # seconds between saves
    delt = delta_t_his/ndiv # time step in seconds

    # this is how we track backwards in time
    if dir_tag == 'reverse':
        delta_t_his = -delta_t_his
        fn_list = fn_list[::-1]
            
    # these lists are used internally to get other variables as needed
    # (they must all be present in the history files)
    vn_list_vel = ['u','v','w']
    vn_list_zh = ['zeta','h']
    vn_list_wind = ['Uwind','Vwind']
    vn_list_other = ['salt', 'temp'] + vn_list_zh + vn_list_vel
    if windage > 0:
        vn_list_other = vn_list_other + vn_list_wind
    # plist_main is what ends up written to output
    plist_main = ['lon', 'lat', 'cs', 'ot', 'z'] + vn_list_other
    
    # Step through times.
    #
    counter_his = 0
    # rot is a list of all the ocean times (sec) in the current file list(e.g. 25)
    # and pot is a single one of these
    for pot in rot[:-1]:

        # get time indices of surrounding (in time) history files
        it0, it1, frt = zfun.get_interpolant(
                np.array(pot), rot, extrap_nan=False)
        # get Datasets
        ds0 = nc4.Dataset(fn_list[it0[0]])
        ds1 = nc4.Dataset(fn_list[it1[0]])
        
        # prepare the fields for nearest neighbor interpolation
        tt0 = time()
        if counter_his == 0:
            h = G['h']
            hf = h[Maskr].data
            if surface == True:
                u0 = ds0['u'][0,-1,:,:]
                u1 = ds1['u'][0,-1,:,:]
                uf0 = u0[Masku].data
                uf1 = u1[Masku].data
                v0 = ds0['v'][0,-1,:,:]
                v1 = ds1['v'][0,-1,:,:]
                vf0 = v0[Maskv].data
                vf1 = v1[Maskv].data
                wf0 = 0; wf1 = 0
                s0 = ds0['salt'][0,-1,:,:]
                s1 = ds1['salt'][0,-1,:,:]
                sf0 = s0[Maskr].data
                sf1 = s1[Maskr].data
                t0 = ds0['temp'][0,-1,:,:]
                t1 = ds1['temp'][0,-1,:,:]
                tf0 = t0[Maskr].data
                tf1 = t1[Maskr].data
                if windage > 0:
                    Uwind0 = ds0['Uwind'][0,:,:]
                    Uwind1 = ds1['Uwind'][0,:,:]
                    Uwindf0 = Uwind0[Maskr].data
                    Uwindf1 = Uwind1[Maskr].data
                    Vwind0 = ds0['Vwind'][0,:,:]
                    Vwind1 = ds1['Vwind'][0,:,:]
                    Vwindf0 = Vwind0[Maskr].data
                    Vwindf1 = Vwind1[Maskr].data
            else:
                u0 = ds0['u'][0,:,:,:]
                u1 = ds1['u'][0,:,:,:]
                uf0 = u0[Masku3].data
                uf1 = u1[Masku3].data
                v0 = ds0['v'][0,:,:,:]
                v1 = ds1['v'][0,:,:,:]
                vf0 = v0[Maskv3].data
                vf1 = v1[Maskv3].data
                w0 = ds0['w'][0,:,:,:]
                w1 = ds1['w'][0,:,:,:]
                wf0 = w0[Maskw3].data
                wf1 = w1[Maskw3].data
                s0 = ds0['salt'][0,:,:,:]
                s1 = ds1['salt'][0,:,:,:]
                sf0 = s0[Maskr3].data
                sf1 = s1[Maskr3].data
                t0 = ds0['temp'][0,:,:,:]
                t1 = ds1['temp'][0,:,:,:]
                tf0 = t0[Maskr3].data
                tf1 = t1[Maskr3].data
                if turb == True:
                    AKs0 = ds0['AKs'][0,:,:,:]
                    AKs1 = ds1['AKs'][0,:,:,:]
                    AKsf0 = AKs0[Maskw3].data
                    AKsf1 = AKs1[Maskw3].data
            z0 = ds0['zeta'][0,:,:]
            z1 = ds1['zeta'][0,:,:]
            zf0 = z0[Maskr].data
            zf1 = z1[Maskr].data
        else:
            # subsequent time steps
            if surface == True:
                u1 = ds1['u'][0,-1,:,:]
                uf0 = uf1.copy()
                uf1 = u1[Masku].data
                v1 = ds1['v'][0,-1,:,:]
                vf0 = vf1.copy()
                vf1 = v1[Maskv].data
                wf0 = 0; wf1 = 0
                s1 = ds1['salt'][0,-1,:,:]
                sf0 = sf1.copy()
                sf1 = s1[Maskr].data
                t1 = ds1['temp'][0,-1,:,:]
                tf0 = tf1.copy()
                tf1 = t1[Maskr].data
                if windage > 0:
                    Uwind1 = ds1['Uwind'][0,:,:]
                    Uwindf0 = Uwindf1
                    Uwindf1 = Uwind1[Maskr].data
                    Vwind1 = ds1['Vwind'][0,:,:]
                    Vwindf0 = Vwindf1
                    Vwindf1 = Vwind1[Maskr].data
            else:
                u1 = ds1['u'][0,:,:,:]
                uf0 = uf1.copy()
                uf1 = u1[Masku3].data
                v1 = ds1['v'][0,:,:,:]
                vf0 = vf1.copy()
                vf1 = v1[Maskv3].data
                w1 = ds1['w'][0,:,:,:]
                wf0 = wf1.copy()
                wf1 = w1[Maskw3].data
                s1 = ds1['salt'][0,:,:,:]
                sf0 = sf1.copy()
                sf1 = s1[Maskr3].data
                t1 = ds1['temp'][0,:,:,:]
                tf0 = tf1.copy()
                tf1 = t1[Maskr3].data
                if turb == True:
                    AKs1 = ds1['AKs'][0,:,:,:]
                    AKsf0 = AKsf1
                    AKsf1 = AKs1[Maskw3].data
            z1 = ds1['zeta'][0,:,:]
            zf0 = zf1.copy()
            zf1 = z1[Maskr].data
        print('Prepare fields for tree %0.4f sec' % (time()-tt0))

        if counter_his == 0:
            if trim_loc == True:
                # remove points on land
                xys = np.array((plon0,plat0)).T
                pmask = maskr[xyT_rho_un.query(xys, n_jobs=-1)[1]]
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
                # NOTE: output is packed in a dict P of arrays ordered as [time, particle]
                P[vn] = np.nan * np.ones((NT,NP))

            # write initial positions to the results arrays
            P['lon'][it0,:] = plon
            P['lat'][it0,:] = plat
            if surface == True:
                pcs[:] = S['Cs_r'][-1]
            P['cs'][it0,:] = pcs
            P['salt'][it0,:] = get_VR_new(sf0, plon, plat, pcs, surface)
            P['temp'][it0,:] = get_VR_new(tf0, plon, plat, pcs, surface)
            V = get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, 0, surface)
            ZH = get_zh_new(zf0,zf1,hf, plon, plat, 0)
            P['u'][it0,:] = V[:,0]
            P['v'][it0,:] = V[:,1]
            P['w'][it0,:] = V[:,2]
            P['zeta'][it0,:] = ZH[:,0]
            P['h'][it0,:] = ZH[:,1]
            P['z'][it0,:] = pcs * ZH.sum(axis=1) + ZH[:,0]
            
        # do the particle tracking for a single pair of history files in ndiv steps
        for nd in range(ndiv):
            
            fr0 = nd/ndiv
            fr1 = (nd + 1)/ndiv
            frmid = (fr0 + fr1)/2
            
            # RK4 integration
            V0 = get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, fr0, surface)
            ZH0 = get_zh_new(zf0,zf1,hf, plon, plat, fr0)
            plon1, plat1, pcs1 = update_position(dxg, dyg, maskr, V0, ZH0, S, delt/2,
                                                 plon, plat, pcs, surface)
                                                 
            V1 = get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon1, plat1, pcs1, frmid, surface)
            ZH1 = get_zh_new(zf0,zf1,hf, plon1, plat1, frmid)
            plon2, plat2, pcs2 = update_position(dxg, dyg, maskr, V1, ZH1, S, delt/2,
                                                 plon, plat, pcs, surface)
                                                 
            V2 = get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon2, plat2, pcs2, frmid, surface)
            ZH2 = get_zh_new(zf0,zf1,hf, plon2, plat2, frmid)
            plon3, plat3, pcs3 = update_position(dxg, dyg, maskr, V2, ZH2, S, delt,
                                                 plon, plat, pcs, surface)
                                                 
            V3 = get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon3, plat3, pcs3, fr1, surface)
            ZH3 = get_zh_new(zf0,zf1,hf, plon3, plat3, fr1)
                                       
            # add windage, calculated from the middle time
            if (surface == True) and (windage > 0):
                Vwind3 = get_wind_new(Uwindf0, Uwindf1, Vwindf0, Vwindf1, plon, plat, frmid, windage)
            else:
                Vwind3 = np.zeros((NP,3))
                
            plon, plat, pcs = update_position(dxg, dyg, maskr, (V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3,
                                              (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                              S, delt, plon, plat, pcs, surface)
            
            # add turbulence to vertical position change (advection already added above)
            if turb == True:
                # pull values of VdAKs and add up to 3-dimensions
                VdAKs = get_dAKs_new(AKsf0, AKsf1, zf0,zf1,hf, plon, plat, pcs, S, frmid)
                
                #VdAKs = get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frmid, surface)
                # print('VdAKs has %d nans out of %d' % (np.isnan(VdAKs).sum(), len(VdAKs)))
                VdAKs3 = np.concatenate((np.zeros((NP,2)), VdAKs[:,np.newaxis]), axis=1)
                
                # update position advecting vertically with 1/2 of AKs gradient
                
                # ZH = get_zh(vn_list_zh, ds0, ds1, plon, plat, pcs, R, frmid, surface)
                # plon_junk, plat_junk, pcs_half = update_position(dxg, dyg, maskr, VdAKs3/2,
                #                 ZH, S, delt, plon, plat, pcs, surface)
                                
                ZH = get_zh_new(zf0,zf1,hf, plon, plat, frmid)
                plon_junk, plat_junk, pcs_half = update_position(dxg, dyg, maskr, VdAKs3/2, ZH, S, delt/2,
                                                     plon, plat, pcs, surface)
                                
                # get AKs at this height, and thence the turbulent velocity
                #Vturb = get_turb(ds0, ds1, VdAKs, delt, plon, plat, pcs_half, R, frmid, surface)
                
                Vturb = get_turb_new(VdAKs, AKsf0, AKsf1, delt, plon, plat, pcs_half, frmid)
                
                # update position vertically for real
                Vturb3 = np.concatenate((np.zeros((NP,2)), Vturb[:,np.newaxis]), axis=1)
                
                # plon_junk, plat_junk, pcs = update_position(dxg, dyg, maskr, Vturb3,
                #                 ZH, S, delt, plon, plat, pcs, surface)
                                
                plon_junk, plat_junk, pcs_half = update_position(dxg, dyg, maskr, Vturb3, ZH, S, delt,
                                                     plon, plat, pcs, surface)
                

        # write positions to the results arrays
        P['lon'][it1,:] = plon
        P['lat'][it1,:] = plat
        if surface == True:
            pcs[:] = S['Cs_r'][-1]
        P['cs'][it1,:] = pcs
        P['salt'][it1,:] = get_VR_new(sf1, plon, plat, pcs, surface)
        P['temp'][it1,:] = get_VR_new(tf1, plon, plat, pcs, surface)
        P['u'][it1,:] = V3[:,0]
        P['v'][it1,:] = V3[:,1]
        P['w'][it1,:] = V3[:,2]
        P['zeta'][it1,:] = ZH3[:,0]
        P['h'][it1,:] = ZH3[:,1]
        P['z'][it1,:] = pcs * ZH3.sum(axis=1) + ZH3[:,0]

        ds0.close()
        ds1.close()

        counter_his += 1

    # by doing this the points are going forward in time
    if dir_tag == 'reverse':
        for vn in plist_main:
            P[vn] = P[vn][::-1,:]

    # and save the time vector (seconds in whatever the model reports)
    P['ot'] = rot

    return P
    
def update_position(dxg, dyg, maskr, V, ZH, S, delta_t, plon, plat, pcs, surface):
    
    # find the new position
    Plon = plon.copy()
    Plat = plat.copy()
    Pcs = pcs.copy()
    # This next step is the actual particle displacement, just done
    # by velocity times a time interval, giving displacements in meters.
    # Each row is a different particle and the columns are (x,y,z) = [0,1,2]
    dX_m = V*delta_t
    # Horizontal advection
    per_m = zfun.earth_rad(Plat)
    clat = np.cos(np.pi*Plat/180.)
    pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
    pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
    Plon += pdx_deg
    Plat += pdy_deg
    
    # Keep particles from being trapped on land
    # minimum grid sizes used when particles approach land boundaries
    xys = np.array((Plon,Plat)).T
    pmask = maskr[xyT_rho_un.query(xys, n_jobs=-1)[1]]
    pcond = pmask < maskr_crit # a Boolean mask
    # move those on land back to their staring point, with a random
    # perturbation of a grid cell to keep them from getting trapped
    if len(pcond) > 0:
        rix = np.random.randint(-1,1,len(plon))
        riy = np.random.randint(-1,1,len(plon))
        Plon[pcond] = plon[pcond] + rix[pcond]*dxg
        Plat[pcond] = plat[pcond] + riy[pcond]*dyg
        # and check again for any stragglers that are still on land
        xys = np.array((Plon,Plat)).T
        pmask = maskr[xyT_rho_un.query(xys, n_jobs=-1)[1]]
        pcond = pmask < maskr_crit
        # setting these back to their starting points (no random perturbation)
        Plon[pcond] = plon[pcond]
        Plat[pcond] = plat[pcond]
        
    # Start of Vertical advection
    H = ZH.sum(axis=1)
    # NOTE: The first column of ZH is the
    # zeta of all particles, and the second is bottom depth h (a positive number),
    # so the H we calculate above is the total water column thickness.
    # The particle z-positions are given by H*pcs + ZH[:,0]
    pdz_s = dX_m[:,2]/H
    # Reflective upper and lower boundary conditions
    if surface == False:
        Pcs_orig = Pcs.copy()
        Pcs += pdz_s
        # Check for bad pcs values.  This was happening a lot becasue of
        # some bugs on the turbulence code that would return bad vertical velocities
        # but since those were fixed this "bad_pcs" diagnostic is all zeros.
        pcs_mask = np.isnan(Pcs)
        if sum(pcs_mask) > 0:
            Pcs[pcs_mask] = Pcs_orig[pcs_mask]
        # Enforce limits on cs.  We use the remainder function to
        # account for cases where the vertical advection may have moved
        # particles more than 2*H 
        hit_top = Pcs > 0
        Pcs[hit_top] = - np.remainder(Pcs[hit_top],1)
        hit_bottom = Pcs < -1
        Pcs[hit_bottom] = -1 - np.remainder(Pcs[hit_bottom],-1)
        # and finally enforce more limits to ensure we always are
        # able to gather salt and other tracer values.
        Pcs[Pcs < S['Cs_r'][0]] = S['Cs_r'][0]
        Pcs[Pcs > S['Cs_r'][-1]] = S['Cs_r'][-1]
        # This is for the cases where particles ended up, for example,
        # in the upper half of the topmost grid cell.  These would
        # be nan'ed out in the get_V call because we specify extrap_nan=True.
    else:
        Pcs[:] = S['Cs_r'][-1]

    return Plon, Plat, Pcs

def get_vel_new(uf0,uf1,vf0,vf1,wf0,wf1, plon, plat, pcs, frac, surface):
    # Get the velocity at all points, at an arbitrary time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1, 0 <= frac <= 1.
    # NOTE: with ndiv=1 this gets called 4 times per hour, or 96 times per day.
    NP = len(plon)
    V = np.zeros((NP,3))
    if surface == True:
        xys = np.array((plon,plat)).T
        # use n_jobs=-1 to use all available cores
        ui0 = uf0[xyT_u.query(xys, n_jobs=-1)[1]]
        vi0 = vf0[xyT_v.query(xys, n_jobs=-1)[1]]
        ui1 = uf1[xyT_u.query(xys, n_jobs=-1)[1]]
        vi1 = vf1[xyT_v.query(xys, n_jobs=-1)[1]]
        ui = (1 - frac)*ui0 + frac*ui1
        vi = (1 - frac)*vi0 + frac*vi1
        V[:,0] = ui
        V[:,1] = vi
    else:
        xys = np.array((plon,plat,pcs)).T
        ui0 = uf0[xyzT_u.query(xys, n_jobs=-1)[1]]
        vi0 = vf0[xyzT_v.query(xys, n_jobs=-1)[1]]
        wi0 = wf0[xyzT_w.query(xys, n_jobs=-1)[1]]
        ui1 = uf1[xyzT_u.query(xys, n_jobs=-1)[1]]
        vi1 = vf1[xyzT_v.query(xys, n_jobs=-1)[1]]
        wi1 = wf1[xyzT_w.query(xys, n_jobs=-1)[1]]
        ui = (1 - frac)*ui0 + frac*ui1
        vi = (1 - frac)*vi0 + frac*vi1
        wi = (1 - frac)*wi0 + frac*wi1
        V[:,0] = ui
        V[:,1] = vi
        V[:,2] = wi
    return V
    
def get_zh_new(zf0,zf1,hf, plon, plat, frac):
    # Get zeta and h at all points, at an arbitrary time between two saves
    NP = len(plon)
    xy = np.array((plon,plat)).T
    zi0 = zf0[xyT_rho.query(xy, n_jobs=-1)[1]]
    zi1 = zf1[xyT_rho.query(xy, n_jobs=-1)[1]]
    hi = hf[xyT_rho.query(xy, n_jobs=-1)[1]]
    zi = (1 - frac)*zi0 + frac*zi1
    ZH = np.zeros((NP,2))
    ZH[:,0] = zi
    ZH[:,1] = hi
    return ZH
    
def get_VR_new(tf, plon, plat, pcs, surface):
    # Get a variable on the z_rho grid at all points.  This only has one time,
    # so it is intended for when we save state variables every hour.
    if surface == True:
        xys = np.array((plon,plat)).T
        ti = tf[xyT_rho.query(xys, n_jobs=-1)[1]]
    else:
        xys = np.array((plon,plat,pcs)).T
        ti = tf[xyzT_rho.query(xys, n_jobs=-1)[1]]
    return ti
    
def get_wind_new(Uwindf0, Uwindf1, Vwindf0, Vwindf1, plon, plat, frac, windage):
    NP = len(plon)
    Vwind3 = np.zeros((NP,3))
    xy = np.array((plon,plat)).T
    Uwind00 = Uwindf0[xyT_rho.query(xy, n_jobs=-1)[1]]
    Uwind11 = Uwindf1[xyT_rho.query(xy, n_jobs=-1)[1]]
    Uwind = (1 - frac)*Uwind00 + frac*Uwind11
    Vwind00 = Vwindf0[xyT_rho.query(xy, n_jobs=-1)[1]]
    Vwind11 = Vwindf1[xyT_rho.query(xy, n_jobs=-1)[1]]
    Vwind = (1 - frac)*Vwind00 + frac*Vwind11
    Vwind3[:,0] = windage*Uwind
    Vwind3[:,1] = windage*Vwind
    return Vwind3
    
def get_AKs_new(AKsf, plon, plat, pcs):
    # Get AKs at all points, at one time.
    xys = np.array((plon,plat,pcs)).T
    AKsi = AKsf[xyzT_w.query(xys, n_jobs=-1)[1]]
    return AKsi
    
def get_dAKs_new(AKsf0, AKsf1, zf0,zf1,hf, plon, plat, pcs, S, frac):
    # create diffusivity gradient for turbulence calculation
    
    # first time
    ZH0 = get_zh_new(zf0,zf1,hf, plon, plat, 0)
    dpcs0 = 1/(ZH0[:,0] + ZH0[:,1]) # change in pcs for a total of a 2m difference
    #     upper variables
    pcs0u = pcs + dpcs0
    pcs0u[pcs0u > S['Cs_w'][-1]] = S['Cs_w'][-1]
    AKs0u = get_AKs_new(AKsf0, plon, plat, pcs0u)
    z0u = pcs0u * ZH0.sum(axis=1)
    #     lower variables
    pcs0b = pcs - dpcs0
    pcs0b[pcs0b < S['Cs_w'][0]] = S['Cs_w'][0]
    AKs0b = get_AKs_new(AKsf0, plon, plat, pcs0b)
    z0b = pcs0b * ZH0.sum(axis=1)
    V0 = (AKs0u-AKs0b)/(z0u-z0b)
    
    # second time
    ZH1 = get_zh_new(zf0,zf1,hf, plon, plat, 1)
    dpcs1 = 1/(ZH1.sum(axis=1)) # change in pcs for a total of a 2m difference
    #     upper variables
    pcs1u = pcs + dpcs1
    pcs1u[pcs1u > S['Cs_w'][-1]] = S['Cs_w'][-1]
    AKs1u = get_AKs_new(AKsf1, plon, plat, pcs1u)
    z1u = pcs1u * ZH1.sum(axis=1)
    #     lower variables
    pcs1b = pcs - dpcs1
    pcs1b[pcs1b < S['Cs_w'][0]] = S['Cs_w'][0]
    AKs1b = get_AKs_new(AKsf1, plon, plat, pcs1b)
    z1b = pcs1b * ZH1.sum(axis=1)
    V1 = (AKs1u-AKs1b)/(z1u-z1b)
    
    # average of times
    V = (1-frac)*V0 + frac*V1
    
    return V

# def get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frac, surface):
#     # create diffusivity gradient for turbulence calculation
#
#     # first time
#     ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R, surface)
#     ZH0[np.isnan(ZH0)] = 0
#     dpcs0 = 1/(ZH0[:,0] + ZH0[:,1]) # change in pcs for a total of a 2m difference
#     #     upper variables
#     pcs0u = pcs + dpcs0
#     pcs0u[pcs0u > S['Cs_w'][-1]] = S['Cs_w'][-1]
#     AKs0u = get_V(['AKs',], ds0, plon, plat, pcs0u, R, surface) # diffusivity
#     z0u = (pcs0u)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
#     #     lower variables
#     pcs0b = pcs - dpcs0
#     pcs0b[pcs0b < S['Cs_w'][0]] = S['Cs_w'][0]
#     AKs0b = get_V(['AKs',], ds0, plon, plat, pcs0b, R, surface) # diffusivity
#     z0b = (pcs0b)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
#     V0 = (AKs0u-AKs0b).squeeze()/(z0u-z0b)
#
#     # second time
#     ZH1 = get_V(vn_list_zh, ds1, plon, plat, pcs, R, surface)
#     ZH1[np.isnan(ZH1)] = 0
#     dpcs1 = 1/(ZH1[:,0] + ZH1[:,1]) # change in pcs for 1m difference
#     #     upper variables
#     pcs1u = pcs + dpcs1
#     pcs1u[pcs1u > S['Cs_w'][-1]] = S['Cs_w'][-1]
#     AKs1u = get_V(['AKs',], ds1, plon, plat, pcs1u, R, surface) # diffusivity
#     z1u = (pcs1u)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
#     #     lower variables
#     pcs1b = pcs - dpcs1
#     pcs1b[pcs1b < S['Cs_w'][0]] = S['Cs_w'][0]
#     AKs1b = get_V(['AKs',], ds1, plon, plat, pcs1b, R, surface) # diffusivity
#     z1b = (pcs1b)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
#     V1 = (AKs1u-AKs1b).squeeze()/(z1u-z1b)
#
#     # average of times
#     V = (1-frac)*V0 + frac*V1
#
#     return V

def get_turb_new(dAKs, AKsf0, AKsf1, delta_t, plon, plat, pcs, frac):
    # get the vertical turbulence correction components
    V0 = get_AKs_new(AKsf0, plon, plat, pcs)
    V1 = get_AKs_new(AKsf1, plon, plat, pcs)
    
    # create weighted average diffusivity
    Vave = (1 - frac)*V0 + frac*V1
    
    # turbulence calculation from Banas, MacCready, and Hickey (2009)
    # w_turbulence = rand*sqrt(2K/dt) + dK/dz
    # rand = random array with normal distribution
    rand = np.random.standard_normal(len(V0))
    V = rand*np.sqrt(2*Vave/delta_t) + dAKs
    
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
    # w_turbulence = rand*sqrt(2K/dt) + dK/dz
    # rand = random array with normal distribution
    rand = np.random.standard_normal(len(V0))
    V = rand*np.sqrt(2*Vave/delta_t) + dAKs
    
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