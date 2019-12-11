"""
This is where you set the run "gtagex" and the initial condition
based on an experiment name passed by the calling code.

"""

import numpy as np

def get_exp_info(exp_name):
    
    EI = {}
    if exp_name == 'fast0':
        # working on tracker_2
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'jdf0'
        
    elif exp_name == 'hc3d':
        # 3D Hood Canal release
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'hc0'
            
    EI['gridname'] = gridname
    EI['tag'] = tag
    EI['ex_name'] = ex_name
    EI['gtagex'] = gridname + '_' + tag + '_' + ex_name
    EI['ic_name'] = ic_name
    return EI
    
def get_ic(ic_name, fn00):
    # routines to set particle initial locations, all numpy arrays
    
    if ic_name == 'hc0': # Hood Canal
        # this fills a volume defined by a rectangular lon, lat region
        # with particles spaced every 10 m from the surface to near the bottom
        lonvec = np.linspace(-123.15, -122.9, 30)
        latvec = np.linspace(47.45, 47.65, 30)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        import zfun
        import zrfun
        G, S, T = zrfun.get_basic_info(fn00)
        hh = zfun.interp2(lonmat, latmat, G['lon_rho'], G['lat_rho'],G['h'])
        plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        hh_vec = hh.flatten()
        for ii in range(len(plon_vec)):
            x = plon_vec[ii]
            y = plat_vec[ii]
            h10 = 10*np.floor(hh_vec[ii]/10) # depth to closest 10 m (above the bottom)
            if h10 >= 10:
                zvec = np.arange(-h10,5,10) # a vector that goes from -h10 to 0 in steps of 10 m
                svec = zvec/hh_vec[ii]
                ns = len(svec)
                if ns > 0:
                    plon00 = np.append(plon00, x*np.ones(ns))
                    plat00 = np.append(plat00, y*np.ones(ns))
                    pcs00 = np.append(pcs00, svec)
    
    else:
        # first create three vectors of initial locations
        # plat00 and plon00 should be the same length,
        # and the length of pcs00 is however many vertical positions you have at
        # each lat, lon (expressed as fraction of depth -1 < pcs < 1)
        if ic_name == 'tn0': # Tacoma Narrows region
            lonvec = np.linspace(-122.65, -122.45, 30)
            latvec = np.linspace(47.2, 47.35, 30)
            lonmat, latmat = np.meshgrid(lonvec, latvec)
            plon_vec = lonmat.flatten()
            plat_vec = latmat.flatten()
            pcs_vec = np.array([-.05])
        
        elif ic_name == 'jdf0': # Mid-Juan de Fuca
            lonvec = np.linspace(-123.85, -123.6, 20)
            latvec = np.linspace(48.2, 48.4, 20)
            lonmat, latmat = np.meshgrid(lonvec, latvec)
            plon_vec = lonmat.flatten()
            plat_vec = latmat.flatten()
            pcs_vec = np.array([-.05])
        # Create full output vectors (each has one value per point).  This
        # code takes each lat, lon location and then assigns it to NSP points
        # corresponding to the vector of pcs values.  However you could write a
        # different version that only released points below a certain depth,
        # or other criterion.
        if len(plon_vec) != len(plat_vec):
            print('WARNING: Problem with length of initial lat, lon vectors')
        NSP = len(pcs_vec)
        NXYP = len(plon_vec)
        plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
        plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
        pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
        plon00 = plon_arr.flatten()
        plat00 = plat_arr.flatten()
        pcs00 = pcs_arr.flatten()
        
    return plon00, plat00, pcs00