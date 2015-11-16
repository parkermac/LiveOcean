"""
Functions for particle tracking.
"""
# setup
import numpy as np
import os; import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path: sys.path.append(alp)
import zfun; reload(zfun) # plotting functions
import netCDF4 as nc4

def get_tracks(fn_list, plon0, plat0, pcs0, delta_t, dir_tag):
 
    plon = plon0.copy()
    plat = plat0.copy()
    pcs = pcs0.copy()
    
    # get basic info
    G, S = zfun.get_basic_info(fn_list[0], getT=False)
       
    # get time vector of history files
    NT = len(fn_list)
    rot = np.nan * np.ones(NT)
    counter = 0
    for fn in fn_list:
        ds = nc4.Dataset(fn)
        rot[counter] = ds.variables['ocean_time'][:].squeeze()
        counter += 1
        ds.close
    
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
    
    # create result arrays
    NP = len(plon)
    Plon = np.nan * np.ones((NT,NP))
    Plat = np.nan * np.ones((NT,NP))
    Pcs = np.nan * np.ones((NT,NP))
          
    Psalt = np.nan * np.ones((NT,NP))
    Ptemp = np.nan * np.ones((NT,NP))
    Pz = np.nan * np.ones((NT,NP))
    Pzeta = np.nan * np.ones((NT,NP))
    Pzbot = np.nan * np.ones((NT,NP))
       
    # Step through times.
    #
    # Currently this is hard-wired to use a 1 hour step with 1 hour saves,
    # first doing a forward step of 1/2 hour, finding the velocity at that time
    # and place, and then using that for a full hour forward step.
    # I think this is what Particulator does, except that it is a longer
    # time step (e.g. 3600 sec instead of 1200 sec).
    #
    counter = 0
    nrot = len(rot)
    for pot in rot[:-1]:
        
        if np.mod(counter,24) == 0:
            print(' - time %d out of %d' % (counter, nrot))
    
        # get time indices  
        iot = zfun.get_interpolant_fast(np.array(pot), rot)
        it0 = iot[0,0].astype(int)
        it1 = iot[0,1].astype(int)
        #print('it0=%d it1=%d' % (it0,it1))
        
        # write positions to the results arrays
        Plon[it0,:] = plon
        Plat[it0,:] = plat
        Pcs[it0,:] = pcs
            
        # get the velocity zeta, and h at all points
        vn_list_vel = ['u','v','w'] 
        vn_list_zh = ['zeta','h']    
        ds0 = nc4.Dataset(fn_list[it0]) 
        V0 = get_V(vn_list_vel, ds0, plon, plat, pcs, R)
        V0[np.isnan(V0)] = 0.0
        ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R)
        #ds0.close()
        
        if counter == 0:
            # get other properties
            vn_list_other = ['salt', 'temp', 'zeta', 'h', ]
            OTH = get_V(vn_list_other, ds0, plon, plat, pcs, R)           
            Psalt[it0,:] = OTH[:,0]
            Ptemp[it0,:] = OTH[:,1]
            this_zeta = OTH[:, vn_list_other.index('zeta')]
            this_h = OTH[:, vn_list_other.index('h')]
            full_depth = this_zeta + this_h
            Pz[it0,:] = pcs * full_depth
            Pzeta[it0,:] = OTH[:,2]
            Pzbot[it0,:] = -OTH[:,3]
       
        # find the new position (half step)
        dX_m = V0*delta_t/2.
        per_m = zfun.earth_rad(plat)
        clat = np.cos(np.pi*plat/180.)
        pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
        pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
        H = ZH0.sum(axis=1)
        pdz_s = dX_m[:,2]/H
        plon_half = plon + pdx_deg
        plat_half = plat + pdy_deg
        pcs_half = pcs + pdz_s
        # enforce limits on cs
        pcs_half[pcs_half < S['Cs_r'][0]] = S['Cs_r'][0]
        pcs_half[pcs_half > S['Cs_r'][-1]] = S['Cs_r'][-1]
        
        # get the velocity zeta, and h at all points
        #ds0 = nc4.Dataset(fn_list[it0]) 
        ds1 = nc4.Dataset(fn_list[it1])
        V0 = get_V(vn_list_vel, ds0, plon_half, plat_half, pcs_half, R)
        V1 = get_V(vn_list_vel, ds1, plon_half, plat_half, pcs_half, R)
        V0[np.isnan(V0)] = 0.0
        V1[np.isnan(V1)] = 0.0
        ZH0 = get_V(vn_list_zh, ds0, plon_half, plat_half, pcs_half, R)
        ZH1 = get_V(vn_list_zh, ds1, plon_half, plat_half, pcs_half, R)
        
        ds0.close()
        #ds1.close()
        
        V = (V0 + V1)/2.
        ZH = (ZH0 + ZH1)/2.
        
        # find the new position
        dX_m = V*delta_t
        per_m = zfun.earth_rad(plat)
        clat = np.cos(np.pi*plat/180.)
        pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
        pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
        H = ZH.sum(axis=1)
        pdz_s = dX_m[:,2]/H
        plon += pdx_deg
        plat += pdy_deg
        pcs += pdz_s
        # enforce limits on cs
        pcs[pcs < S['Cs_r'][0]] = S['Cs_r'][0]
        pcs[pcs > S['Cs_r'][-1]] = S['Cs_r'][-1]
        
        # find properties at the new position
        OTH = get_V(vn_list_other, ds1, plon, plat, pcs, R)
        Psalt[it1,:] = OTH[:,0]
        Ptemp[it1,:] = OTH[:,1]
        this_zeta = OTH[:, vn_list_other.index('zeta')]
        this_h = OTH[:, vn_list_other.index('h')]
        full_depth = this_zeta + this_h
        Pz[it1,:] = pcs * full_depth
        Pzeta[it1,:] = OTH[:,2]
        Pzbot[it1,:] = -OTH[:,3]
        
        ds1.close()
        
        counter += 1
    
    # add final point    
    Plon[it1,:] = plon
    Plat[it1,:] = plat
    Pcs[it1,:] = pcs
    
    # by doing this the points are going forward in time
    if dir_tag == 'reverse':
        Plon = Plon[::-1,:]
        Plat = Plat[::-1,:]
        Pcs = Pcs[::-1,:]
        
        Psalt = Psalt[::-1,:]
        Ptemp = Ptemp[::-1,:]
        Pz = Pz[::-1,:]
        Pzeta = Pzeta[::-1,:]
        Pzbot = Pzbot[::-1,:]
            
    # and the time vector (seconds in whatever the model reports)
    Pot = rot
    
    # put it all in a dict
    P = dict()
    P['lon'] = Plon
    P['lat'] = Plat
    P['cs'] = Pcs
    P['ot'] = Pot
    P['salt'] = Psalt
    P['temp'] = Ptemp
    P['z'] = Pz
    P['zeta'] = Pzeta
    P['zbot'] = Pzbot
    
    return P, G, S

def fix_masked(V,mask_val):    
    try:
        if V.mask.any():
            if mask_val == 99:
                vm = np.array(~V.mask, dtype=bool).flatten()
                vf = V.flatten()
                vgood = vf[vm].mean()
                V[V.mask] = vgood
            else:
                V[V.mask] = mask_val
            V = V.data
    except AttributeError:
            pass
    return V
    
def get_V(vn_list, ds, plon, plat, pcs, R):

    # get interpolant arrays
    ilonr_a = zfun.get_interpolant_fast(plon, R['rlonr'])
    ilatr_a = zfun.get_interpolant_fast(plat, R['rlatr'])
    ilonu_a = zfun.get_interpolant_fast(plon, R['rlonu'])
    ilatu_a = zfun.get_interpolant_fast(plat, R['rlatu'])
    ilonv_a = zfun.get_interpolant_fast(plon, R['rlonv'])
    ilatv_a = zfun.get_interpolant_fast(plat, R['rlatv'])
    icsr_a = zfun.get_interpolant_fast(pcs, R['rcsr'])
    icsw_a = zfun.get_interpolant_fast(pcs, R['rcsw'])
        
    NV = len(vn_list)
    NP = len(plon)
    
    # get a interpolated values    
    V = np.nan * np.ones((NP,NV))
    for ip in range(NP):
        vcount = 0
        for vn in vn_list:    
            if vn in ['salt','temp','zeta','h']:
                ics = icsr_a[ip,:]
                ilat = ilatr_a[ip,:]
                ilon = ilonr_a[ip,:]
            elif vn in ['u']:
                ics = icsr_a[ip,:]
                ilat = ilatu_a[ip,:]
                ilon = ilonu_a[ip,:]
            elif vn in ['v']:
                ics = icsr_a[ip,:]
                ilat = ilatv_a[ip,:]
                ilon = ilonv_a[ip,:]
            elif vn in ['w']:
                ics = icsw_a[ip,:]
                ilat = ilatr_a[ip,:]
                ilon = ilonr_a[ip,:]
            if vn in ['salt','temp','u','v','w']:      
                box0 = ds.variables[vn][0,ics[:2].astype(int),
                    ilat[:2].astype(int),ilon[:2].astype(int)].squeeze()
                if vn in ['u', 'v', 'w']:
                    box0 = fix_masked(box0, 0.0)
                elif vn in ['salt','temp']:
                    #print('box0')
                    #print(box0)
                    box0 = fix_masked(box0, 99)
                    #print('box0 fixed')
                    #print(box0)
                az = np.array([1-ics[2], ics[2]])        
                ay = np.array([1-ilat[2], ilat[2]]).reshape((1,2))
                ax = np.array([1-ilon[2], ilon[2]]).reshape((1,1,2))                            
                V[ip,vcount] = (az*( ( ay*((ax*box0).sum(-1)) ).sum(-1) )).sum()
            elif vn in ['zeta']:      
                box0 = ds.variables[vn][0,
                    ilat[:2].astype(int),ilon[:2].astype(int)].squeeze()
                az = 1.
                ay = np.array([1-ilat[2], ilat[2]]).reshape((1,2))
                ax = np.array([1-ilon[2], ilon[2]]).reshape((1,1,2))                            
                V[ip,vcount] = (az*( ( ay*((ax*box0).sum(-1)) ).sum(-1) )).sum()
            elif vn in ['h']:      
                box0 = ds.variables[vn][
                    ilat[:2].astype(int),ilon[:2].astype(int)].squeeze()
                az = 1.
                ay = np.array([1-ilat[2], ilat[2]]).reshape((1,2))
                ax = np.array([1-ilon[2], ilon[2]]).reshape((1,1,2))                            
                V[ip,vcount] = (az*( ( ay*((ax*box0).sum(-1)) ).sum(-1) )).sum()
            vcount += 1
    return V