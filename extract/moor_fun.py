"""
Mooring extraction functions.
"""

import numpy as np
import netCDF4 as nc
import zfun

def get_sta_dict(job_name):
    if job_name == 'kd_array': # Kathy Donohue, Slow Slip related 2018_09
        sta_dict = {
            'kda00': (-125.155750, 44.447250),
            'kda01': (-125.206174, 44.447359),
            'kda02': (-125.256597, 44.447469),
            'kda03': (-125.307021, 44.447578),
            'kda04': (-125.357446, 44.447688),
            'kda05': (-125.407870, 44.447797),
            'kda06': (-125.458295, 44.447906),
            'kda07': (-125.155750, 44.402250),
            'kda08': (-125.206134, 44.401207),
            'kda09': (-125.256516, 44.400165),
            'kda10': (-125.306896, 44.399123),
            'kda11': (-125.357275, 44.398081),
            'kda12': (-125.407651, 44.397039),
            'kda13': (-125.458026, 44.395997),
            'kda14': (-125.155750, 44.357249),
            'kda15': (-125.206095, 44.355605),
            'kda16': (-125.256436, 44.353961),
            'kda17': (-125.306775, 44.352317),
            'kda18': (-125.357111, 44.350672),
            'kda19': (-125.407445, 44.349028),
            'kda20': (-125.457775, 44.347384),
            'kda21': (-125.085750, 44.177247),
            'kda22': (-125.135941, 44.175703),
            'kda23': (-125.186129, 44.174159),
            'kda24': (-125.236314, 44.172614),
            'kda25': (-125.286497, 44.171070),
            'kda26': (-125.336677, 44.169525),
            'kda27': (-125.386855, 44.167981)
            }
        v2_list = ['zeta', 'Pair']
        v3_list_rho = ['rho', 'u', 'v']
        v3_list_w = []
    return sta_dict, v2_list, v3_list_rho, v3_list_w

def get_itp_dict(sta_dict, G):
    # get interpolants
    itp_dict = dict()
    for sta_name in sta_dict.keys():
        # target position
        Lon = np.array(sta_dict[sta_name][0])
        Lat = np.array(sta_dict[sta_name][1])
        # get interpolants for this point
        Xi0 = dict(); Yi0 = dict()
        Xi1 = dict(); Yi1 = dict()
        Aix = dict(); Aiy = dict()
        for grd in ['rho', 'u', 'v']:
            xx = G['lon_' + grd][1,:]
            yy = G['lat_' + grd][:,1]
            xi0, xi1, xfr = zfun.get_interpolant(Lon, xx, extrap_nan=True)
            yi0, yi1, yfr = zfun.get_interpolant(Lat, yy, extrap_nan=True)
            Xi0[grd] = xi0
            Yi0[grd] = yi0
            Xi1[grd] = xi1
            Yi1[grd] = yi1
            # create little arrays that are used in the actual interpolation
            Aix[grd] = np.array([1-xfr, xfr]).reshape((1,1,2))
            Aiy[grd] = np.array([1-yfr, yfr]).reshape((1,2))
        itp_dict[sta_name] = (Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
        
    return itp_dict

def get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy):
    dims = ds.variables[vv].dimensions
    if 'eta_rho' in dims:
        grd = 'rho'
    elif 'eta_u' in dims:
        grd = 'u'
    elif 'eta_v' in dims:
        grd = 'v'
    else:
        print('grid error!')
    xi0 = Xi0[grd]; yi0 = Yi0[grd]
    xi1 = Xi1[grd]; yi1 = Yi1[grd]
    aix = Aix[grd]; aiy = Aiy[grd]

    xi01 = np.array([xi0, xi1]).flatten()
    yi01 = np.array([yi0, yi1]).flatten()
    return xi01, yi01, aix, aiy
    
def start_netcdf(out_fn, N, NT, v0_list, v1_list, v2_list, v3_list_rho, v3_list_w, V_long_name, V_units):

    foo = nc.Dataset(out_fn, 'w')
    foo.createDimension('scalar', 1)
    foo.createDimension('s_rho', N)
    foo.createDimension('s_w', N+1)
    foo.createDimension('ocean_time', NT)

    for vv in v0_list:
        v_var = foo.createVariable(vv, float, ('scalar'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v1_list:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v2_list:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v3_list_rho:
        v_var = foo.createVariable(vv, float, ('ocean_time','s_rho'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]
    for vv in v3_list_w:
        v_var = foo.createVariable(vv, float, ('ocean_time', 's_w'))
        v_var.long_name = V_long_name[vv]
        v_var.units = V_units[vv]

    foo.close()