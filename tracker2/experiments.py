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
        
    elif exp_name == 'fast1':
        # working on tracker_2
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'tn0'
    
    elif exp_name == 'hc3d':
        # 3D Hood Canal release
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'hc1'
        
    elif exp_name == 'trap0':
        # test of particle trapping in rivers
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'skok'
        
    elif exp_name == 'ai0':
        # test of particle trapping in rivers
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'sea'

    elif exp_name == 'vmix':
        # test of vertical mixing; use with the new flag:
        # -no_advection True, so a full command would be
        # python tracker.py -exp vmix -3d True -clb True -no_advection True
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
        ic_name = 'vm1'
            
            
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
                    
    elif ic_name == 'hc1': # Hood Canal using TEF "segments"
        # NOTE: this one requires information about grid indices in certain regions of
        # the Salish Sea, and is very specific to the analysis in x_tef.
        # Not for general use.
        import pickle
        # select the indir
        import Lfun
        Ldir = Lfun.Lstart()
        indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
        indir = indir0 + 'flux/'
        # load data
        ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))

        import zrfun
        G = zrfun.get_basic_info(fn00, only_G=True)
        h = G['h']
        xp = G['lon_rho']
        yp = G['lat_rho']
        
        plon_vec = np.array([])
        plat_vec = np.array([])
        hh_vec = np.array([])
        
        seg_list = ['H3','H4','H5','H6','H7','H8']
        for seg_name in seg_list:
            ji_seg = ji_dict[seg_name]
            for ji in ji_seg:
                plon_vec = np.append(plon_vec, xp[ji])
                plat_vec = np.append(plat_vec, yp[ji])
                hh_vec = np.append(hh_vec, h[ji])
                    
        plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
    
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
                    
        # trim the length of the I.C. vectors to a limited number
        NP = len(plon00)
        npmax = 10000
        nstep = max(1,int(NP/npmax))
        plon00 = plon00[::nstep]
        plat00 = plat00[::nstep]
        pcs00 = pcs00[::nstep]
        
    elif ic_name == 'sea': # mouth of AI
        import pickle
        # get the ji_dict which has indices of named segments
        import Lfun
        Ldir = Lfun.Lstart()
        ji_dict = pickle.load(open(Ldir['LOo'] + 'tef/volumes_cas6/ji_dict.p', 'rb'))

        import zrfun
        G = zrfun.get_basic_info(fn00, only_G=True)
        xp = G['lon_rho']; yp = G['lat_rho']; h = G['h']
        plon_vec = np.array([]); plat_vec = np.array([]); hh_vec = np.array([])
        for seg_name in ['J4']:
            ji_seg = ji_dict[seg_name]
            for ji in ji_seg: # ji is a tuple
                plon_vec = np.append(plon_vec, xp[ji])
                plat_vec = np.append(plat_vec, yp[ji])
                hh_vec = np.append(hh_vec, h[ji])
                    
        plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
        for ii in range(len(plon_vec)):
            x = plon_vec[ii]
            y = plat_vec[ii]
            DZ = 5
            hdz = DZ*np.floor(hh_vec[ii]/DZ) # depth to closest DZ m (above the bottom)
            if hdz >= DZ:
                zvec = np.arange(-hdz,DZ,DZ) # a vector that goes from -hdz to 0 in steps of DZ m
                svec = zvec/hh_vec[ii]
                ns = len(svec)
                if ns > 0:
                    plon00 = np.append(plon00, x*np.ones(ns))
                    plat00 = np.append(plat00, y*np.ones(ns))
                    pcs00 = np.append(pcs00, svec)
        
        # trim the length of the I.C. vectors to a limited number
        NP = len(plon00)
        print('-- NP = %d (original)' % (NP))
        npmax = 10000
        if NP > npmax:
            nstep = max(1,int(np.floor(NP/npmax)))
            plon00 = plon00[::nstep]
            plat00 = plat00[::nstep]
            pcs00 = pcs00[::nstep]
            print('-- NP = %d (trimmed)' % (len(pcs00)))
        else:
            print('-- NP = %d' % (len(pcs00)))

                    
    # cases that use ic_from_meshgrid
    elif ic_name == 'tn0': # Tacoma Narrows region
        lonvec = np.linspace(-122.65, -122.45, 30)
        latvec = np.linspace(47.2, 47.35, 30)
        pcs_vec = np.array([-.05])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
    elif ic_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        pcs_vec = np.array([-.05])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
    elif ic_name == 'skok': # head of the Skokomish River
        lonvec = np.linspace(-123.171, -123.163, 10)
        latvec = np.linspace(47.306, 47.312, 10)
        pcs_vec = np.array([-.05])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
    elif ic_name == 'vm1': # for test of vertical mixing
        lonvec = np.array([-125.35, -124.0, -122.581]); latvec = np.array([47.847, 48.3, 48.244])
        # slope off JdF, middle of JdF, Whidbey Basin
        import zrfun
        S = zrfun.get_basic_info(fn00, only_S=True)
        print(S['Cs_r'][0])
        print(S['Cs_r'][-1])
        pcs_vec = np.linspace(S['Cs_r'][0],S['Cs_r'][-1],num=1000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
        
    return plon00, plat00, pcs00
    
def ic_from_meshgrid(lonvec, latvec, pcs_vec):
    # First create three vectors of initial locations (as done in some cases above).
    # plat00 and plon00 should be the same length, and the length of pcs00 is
    # however many vertical positions you have at each lat, lon
    # (expressed as fraction of depth -1 < pcs < 1)
    #
    # Then here we create full output vectors (each has one value per point).
    # This code takes each lat, lon location and then assigns it to NSP points
    # corresponding to the vector of pcs values.
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
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
    
def ic_from_list(lonvec, latvec, pcs_vec):
    # Like ic_from_meshgrid() but treats the lon, lat lists like lists of mooring locations.
    plon_vec = lonvec
    plat_vec = latvec
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00
    