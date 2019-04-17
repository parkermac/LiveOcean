"""
Functions for the carbon code.
"""
import scipy.io as io
import netCDF4 as nc
import subprocess
import numpy as np
import time
import seawater as sw
import sys

# assume path to alpha set by calling function
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun

def fillit(a):
    # ensures a is an array with nan's for masked values
    # instead of a masked array
    if isinstance(a, np.ma.MaskedArray):
        a = a.filled(np.nan)
    return a
    
def get_full(fn):
    # gets v_dict for the full extent of one history file
    ds = nc.Dataset(fn)
    G, S, T = zrfun.get_basic_info(fn)        

    # extract needed info from history file
    v_dict = dict()
    if print_info:
        print('\nINPUT Variable Info:')
    for vn in ['alkalinity', 'TIC', 'salt', 'temp','rho']:
        v = ds[vn][:]
        v = fillit(v)
        v_dict[vn] = v
                
    # create depth, pressure, and in situ temperature
    h = ds['h'][:]
    h = fillit(h)
    lat = G['lat_rho'][:]
    z_rho = zrfun.get_z(h, 0*h, S, only_rho=True)
    depth = -z_rho
    pres = sw.pres(depth, lat)
    v_dict['pres'] = pres
    temp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, v_dict['pres'])
    v_dict['temp'] = temp
    
    # convert from umol/L to umol/kg using in situ dentity
    v_dict['alkalinity'] = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
    v_dict['TIC'] = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
    
    # clean up
    v_dict.pop('rho') # no longer needed, so don't pass to worker
    ds.close()
    
    return v_dict

def get_multi_layers(fn_list, NZ, aa):
    ii = 0
    tt0 = time.time()
    for fn in fn_list:
        # get the fields
        v_dict, plon, plat = get_layer(fn, NZ=NZ, aa=aa)
        if ii == 0:
            vv_dict = v_dict.copy()
            for vn in v_dict.keys():
                v = v_dict[vn]
                NR, NC = v.shape
                vv_dict[vn] = v.reshape((1, NR, NC))
        else:
            for vn in v_dict.keys():
                v = v_dict[vn]
                NR, NC = v.shape
                vv_dict[vn] = np.concatenate((vv_dict[vn], v.reshape((1, NR, NC))), axis=0)
        ii += 1
    print('\nTIME to gather and process input fields:')
    print(' -- %d input files' % (len(fn_list)))
    print(' -- %0.1f seconds' % (time.time() - tt0))
    sys.stdout.flush()
    
    return vv_dict, plon, plat

def get_layer(fn, NZ=-1, aa=[], print_info=False):
    # function to extract and process fields from a history file
    # returning a dict of arrays that can be passed to CO2SYS.m
    
    # default is to get full surface field
    
    ds = nc.Dataset(fn)
    G, S, T = zrfun.get_basic_info(fn)
    if len(aa) == 4:
        # find indices that encompass region aa
        i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
        i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
        j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
        j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
    else:
        # full region
        i0 = 0; j0 = 0
        j1, i1 = G['lon_rho'].shape
        i1 += 1; j1 += 1
        
    plon = G['lon_psi'][j0:j1-1, i0:i1-1]
    plat = G['lat_psi'][j0:j1-1, i0:i1-1]

    # extract needed info from history file
    v_dict = dict()
    if print_info:
        print('\nINPUT Variable Info:')
    for vn in ['alkalinity', 'TIC', 'salt', 'temp','rho']:
        v = ds[vn][0,NZ, j0:j1, i0:i1]
        v = fillit(v)
        v_dict[vn] = v
        if print_info:
            name = ds[vn].long_name
            try:
                units = ds[vn].units
            except AttributeError:
                units = ''
            vmax = np.nanmax(v)
            vmin = np.nanmin(v)
            v_dict[vn] = v
            print('%25s (%25s) max = %6.1f min = %6.1f' %
                (name, units, vmax, vmin))
                
    # create depth, pressure, and in situ temperature
    h = ds['h'][j0:j1, i0:i1]
    h = fillit(h)
    lat = G['lat_rho'][j0:j1, i0:i1]
    z_rho = zrfun.get_z(h, 0*h, S, only_rho=True)
    depth = -z_rho[NZ, :, :].squeeze()
    pres = sw.pres(depth, lat)
    v_dict['pres'] = pres
    # assume potential temperature is close enough
    # temp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, v_dict['pres'])
    # v_dict['temp'] = temp
    
    # convert from umol/L to umol/kg using in situ dentity
    v_dict['alkalinity'] = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
    v_dict['TIC'] = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
    
    # clean up
    v_dict.pop('rho') # no longer needed, so don't pass to worker
    ds.close()
    
    return v_dict, plon, plat
    
def make_carbon(v_dict, tempdir, print_info=False):
    # pass data to CO2SYS.m

    # save v_dict to temporary file for matlab
    Lfun.make_dir(tempdir, clean=True)
    in_fn = tempdir + 'worker_input.mat'
    io.savemat(in_fn,v_dict)
    
    # run the matlab code
    if print_info:
        tt0 = time.time()
    func = "worker(\'" + tempdir + "\')"
    cmd = Ldir['which_matlab']
    run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
    proc = subprocess.run(run_cmd, stdout=subprocess.PIPE, cwd=Ldir['LO']+'co2sys')
    if print_info:
        print('\n TIME to do calculation:')
        print(' -- run CO2SYS %0.1f s' % (time.time()-tt0))

def get_carbon(tempdir):
    # load the output of CO2SYS
    PH = io.loadmat(tempdir + 'PH.mat')['PH']
    OM = io.loadmat(tempdir + 'OM.mat')['OM']
    
    return PH, OM