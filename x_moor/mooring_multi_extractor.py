"""
Extract multiple mooring-like records.
The input is structured to conform with layer_extractor.py.

"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun
import netCDF4 as nc

import moor_fun as mfun
from importlib import reload
reload(mfun)

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.07.20')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.07.22')
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default='hourly')
args = parser.parse_args()

# get list of history files to plot
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
fn_list = Lfun.get_fn_list(args.list_type, Ldir, args.date_string0, args.date_string1)

# specify list of moorings to work on
if False:
    # specify by hand
    sta_dict = {'m1': (-124.6, 48.46),
                'm2': (-124.6, 47)}
    limit_lists = False
else:
    # read in a custom job
    job_name = 'willapa_bc'
    sta_dict, v2_list, v3_list_rho, v3_list_w = mfun.get_sta_dict(job_name)
    limit_lists = False
    
    # job_name = 'kd_array'
    # sta_dict, v2_list, v3_list_rho, v3_list_w = mfun.get_sta_dict(job_name)
    # limit_lists = True

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

# get interpolants
itp_dict = mfun.get_itp_dict(sta_dict, G)

# INITIALIZATION
ds = nc.Dataset(fn)
# generating some lists
v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']
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
            and (vv != 'ocean_time') ):
            v2_list.append(vv)
        elif ( ('ocean_time' in vdim) and ('s_rho' in vdim) ):
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
# save static variables
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    Xi0, Yi0, Xi1, Yi1, Aix, Aiy = itp_dict[sta_name]
    for vv in v0_list:
        xi01, yi01, aix, aiy = mfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
        vvtemp = ds.variables[vv][yi01, xi01].squeeze()
        foo = nc.Dataset(out_fn, 'a')
        foo[vv][:] =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
        foo.close()
ds.close()
# ENF OF INITIALIZATION

# EXTRACT TIME-DEPENDENT FIELDS
count = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    
    if np.mod(count,24)==0:
        print(' working on %d of %d' % (count, NT))
        sys.stdout.flush()
        
    for sta_name in sta_dict.keys():
        out_fn = out_fn_dict[sta_name]
        foo = nc.Dataset(out_fn, 'a')
        Xi0, Yi0, Xi1, Yi1, Aix, Aiy = itp_dict[sta_name]
        for vv in v1_list:
            vtemp = ds.variables[vv][:].squeeze()
            foo[vv][count] = vtemp
        for vv in v2_list:
            xi01, yi01, aix, aiy = mfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count] = vtemp
        for vv in v3_list_rho:
            xi01, yi01, aix, aiy = mfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count,:] = vtemp
        for vv in v3_list_w:
            xi01, yi01, aix, aiy = mfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count,:] = vtemp
        foo.close()
    count += 1

    ds.close()
# END OF EXTRACTING TIME-DEPENDENT FIELDS

# create z_rho and z_w (has to be done after we have zeta)
for sta_name in sta_dict.keys():
    out_fn = out_fn_dict[sta_name]
    foo = nc.Dataset(out_fn, 'a')

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
    
    foo.close()

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


