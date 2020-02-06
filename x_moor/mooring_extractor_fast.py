"""
Extract multiple mooring-like records.
The input is structured to conform with layer_extractor.py.

Should be a bit faster than mooring_extractor.py because we
skip all the interpolation and just use a model cell that is
near the station, and don't bother interpolating u, and v perfectly.

NOTE: when specifying a negative longitude or latitude from the command line
you have to do it like:
        -lon " -124.5"
including the space before the minus sign.  Otherwise argparse gets confused.

Example call from the command line:

run mooring_extractor.py -sn Dabob -lon " -122.85" -lat 47.7

"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun
import netCDF4 as nc

from time import time

from importlib import reload

import moor_fun as mfun
reload(mfun)

# Load the module that defines mooring lists.  First look for a
# "user_" one that you can put in this directory, copying the
# format in moor_lists.py.
if os.path.isfile('user_moor_lists.py'):
    import user_moor_lists as ml
else:
    import moor_lists as ml
reload(ml)

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2019.07.04')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2019.07.05')
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default='hourly')
parser.add_argument('-v', '--verbose', default=False, type=boolean_string)
# see alpha/Lfun.get_fn_list() for acceptable list_type values

# Mooring arguments.  You MUST supply all arguments for either (1) or (2)
#
# (1) Define a single mooring from the command line.
parser.add_argument('-sn', '--sta_name', nargs='?', type=str, default='blank')
parser.add_argument('-lon', '--lon_str', nargs='?', type=str, default='')
parser.add_argument('-lat', '--lat_str', nargs='?', type=str, default='')
#
# (2) Specify the job name to use in mfun.get_sta_dict(job_name).
parser.add_argument('-jn', '--job_name', nargs='?', type=str, default='blank')

args = parser.parse_args()
verbose = args.verbose

# get list of history files to plot
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
fn_list = Lfun.get_fn_list(args.list_type, Ldir, args.date_string0, args.date_string1)

# specify list of moorings to work on

limit_lists = False # default is to get all variables
# Initialize lists of varibles to get.  If you leave them empty
# then the calling code will get everything.
v2_list = [] # 2-D variables (like 'zeta')
v3_list_rho = [] # 3-D variables on the vertical rho grid (like 'salt')
v3_list_w = [] # 3-D varibles on the vertical w grid (like 'w')

if (args.sta_name != 'blank') and (args.lon_str != '') and (args.lat_str != ''):
    # First see if we should make info for single station
    # from command line arguments.
    sta_dict = {args.sta_name: (float(args.lon_str), float(args.lat_str))}
elif args.job_name != 'blank':
    # If not then assume we will
    # use a job defined in moor_lists.py (or user_moor_lists.py if it exists)
    sta_dict, v2_list, v3_list_rho, v3_list_w = ml.get_sta_dict(args.job_name)
else:
    print('Inconsistent command line arguments for mooring.')
    sys.exit()
    
if (len(v2_list)>0) or (len(v3_list_rho)>0) or (len(v3_list_w)>0):
    limit_lists = True
    
# name output files
out_fn_dict = dict()
for sta_name in sta_dict.keys():
    # make sure the output directory exists
    outdir0 = Ldir['LOo'] + 'moor/'
    Lfun.make_dir(outdir0)
    mod_string = (Ldir['gtagex'] + '_' + args.date_string0 + '_' + args.date_string1)
    outdir = outdir0 + mod_string + '/'
    Lfun.make_dir(outdir)
    # name output file
    out_fn = (outdir + sta_name + '_' +  args.list_type + '.nc')
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass
    out_fn_dict[sta_name] = out_fn

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
N = S['N']
NT = len(fn_list)

lon = G['lon_rho']
lat = G['lat_rho']
mask = G['mask_rho']
xvec = lon[0,:].flatten()
yvec = lat[:,0].flatten()

# automatically correct for moorings that are on land
# Note that this only checks rho points - what if a u or v point is masked?
ji_dict = {}
for sta in sta_dict.keys():
    xy = sta_dict[sta]
    slon = xy[0]
    slat = xy[1]
    i0 = zfun.find_nearest_ind(xvec, slon)
    j0 = zfun.find_nearest_ind(yvec, slat)
    # find indices of nearest good point
    # Note: tha mask in G is Boolean, not 0, 1, but testing with 1 or 0 works.
    if mask[j0,i0] == 1:
        print(sta + ': point ok')
    elif mask[j0,i0] == 0:
        print(sta + ':point masked')
        i0, j0 = mfun.get_ij_good(lon, lat, mask, xvec, yvec, i0, j0)
        newlon = xvec[i0]
        newlat = yvec[j0]
        print(' Replacing (%0.3f,%0.3f) with (%0.3f,%0.3f)' % (slon, slat, newlon, newlat))
        sta_dict[sta] = (newlon, newlat)
    ji_dict[sta] = (j0, i0)

tt0 = time()
# INITIALIZATION
ds = nc.Dataset(fn,'r')
# generating some lists
v0_list = ['h', 'lon_rho', 'lat_rho']
v1_list = ['ocean_time']
if limit_lists == True:
    pass # use definitions from above
else:
    v2_list = []
    v3_list_rho = []
    v3_list_w = []
    for vv in ds.variables:
        vdim = ds.variables[vv].dimensions
        if ( ('ocean_time' in vdim)
            and ('s_rho' not in vdim)
            and ('s_w' not in vdim)
            and (vv != 'ocean_time') and (vv != 'ssflux')):
            v2_list.append(vv)
        elif ( ('ocean_time' in vdim) and ('s_rho' in vdim) and (vv != 'CaCO3')):
            v3_list_rho.append(vv)
        elif ( ('ocean_time' in vdim) and ('s_w' in vdim) ):
            v3_list_w.append(vv)

# generate dictionaries of long names and units
V_long_name = dict()
V_units = dict()
v_all_list = v0_list + v1_list + v2_list + v3_list_rho + v3_list_w
for vv in v_all_list:
    try:
        V_long_name[vv] = ds.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = ds.variables[vv].units
    except:
        V_units[vv] = ''
# start netcdf files
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    mfun.start_netcdf(out_fn, N, NT, v0_list, v1_list, v2_list,
        v3_list_rho, v3_list_w, V_long_name, V_units)
        
# open all output files into a dict
foo_dict = {}
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    foo = nc.Dataset(out_fn, 'a')
    foo_dict[sta_name] = foo
        
# save static variables
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    for vv in v0_list:
        j0, i0 = ji_dict[sta_name]
        foo = foo_dict[sta_name]
        foo[vv][:] = ds[vv][j0, i0]
ds.close()

# END OF INITIALIZATION
if verbose:
    print(' -- initialization took %0.2f seconds' % (time()-tt0))

# EXTRACT TIME-DEPENDENT FIELDS


count = 0
for fn in fn_list:
    
    tt1 = time()
    
    ds = nc.Dataset(fn, 'r')
    
    if np.mod(count,24)==0:
        print(' working on %d of %d' % (count, NT))
        sys.stdout.flush()

    for sta_name in sta_dict.keys():
        j0, i0 = ji_dict[sta_name]
        tt2 = time()
        foo = foo_dict[sta_name]
        for vv in v1_list:
            val = ds[vv][0]
            foo[vv][count] = val.data
        for vv in v2_list:
            val = ds[vv][0, j0, i0]
            foo[vv][count] = val.data
        # for vv in v3_list_rho + v3_list_w:
        #     val = ds[vv][0,:,j0,i0]
        #     #foo[vv][count,:] = val
        if verbose:
            print(' ... station took %0.2f seconds' % (time()-tt2))
    count += 1
    ds.close()
    
    if verbose:
        print(' -- this history file took %0.2f seconds' % (time()-tt1))
    

# END OF EXTRACTING TIME-DEPENDENT FIELDS

# create z_rho and z_w (has to be done after we have zeta)
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    
    foo = foo_dict[sta_name]
    
    #foo = nc.Dataset(out_fn, 'a')

    zeta = foo['zeta'][:].squeeze()
    hh = foo['h'][:] * np.ones_like(zeta)
    z_rho, z_w = zrfun.get_z(hh, zeta, S)
    
    v_var = foo.createVariable('z_rho', float, ('ocean_time','s_rho'))
    v_var.long_name = 'z on rho points (positive up)'
    v_var.units = 'm'
    v_var[:] = z_rho.T
    
    v_var = foo.createVariable('z_w', float, ('ocean_time','s_w'))
    v_var.long_name = 'z on w points (positive up)'
    v_var.units = 'm'
    v_var[:] = z_w.T
    
    #foo.close()
# close all output files
for sta_name in sta_dict.keys():
    foo_dict[sta_name].close()

# finale
import collections
result_dict = collections.OrderedDict()
result_dict['out_fn'] = out_fn
result_dict['date_string0'] = args.date_string0
result_dict['date_string1'] = args.date_string1
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'
print('')
for k in result_dict.keys():
    print('%s: %s' % (k, result_dict[k]))


