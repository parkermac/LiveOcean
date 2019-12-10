"""
This is where you set the run "gtagex" and the initial condition
based on an experiment name passed by the calling code.

"""

import numpy as np

def get_exp_info(exp_name):
    
    EI = {}
    if exp_name == 'fast0': # working on tracker_2
            gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
            ic_name = 'jdf0'
            
    EI['gridname'] = gridname
    EI['tag'] = tag
    EI['ex_name'] = ex_name
    EI['gtagex'] = gridname + '_' + tag + '_' + ex_name
    EI['ic_name'] = ic_name
    return EI
    
def get_ic(ic_name, fn00):
    
    # routines to set particle initial locations, all numpy arrays
    #
    # first create three vectors of initial locations
    # plat00 and plon00 should be the same length,
    # and the length of pcs00 is however many vertical positions you have at
    # each lat, lon (expressed as fraction of depth -1 < pcs < 1)
    
    if ic_name == 'hc0': # Hood Canal
        lonvec = np.linspace(-122.65, -122.45, 30)
        latvec = np.linspace(47.2, 47.35, 30)
        lonmat, latmat = np.meshgrid(lonvec, latvec)
        plon_vec = lonmat.flatten()
        plat_vec = latmat.flatten()
        pcs_vec = np.array([-.05])

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