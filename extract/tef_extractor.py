"""
Extract fields at a section which may be used later for TEF analysis
of transport and transport-weighted properties.

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

# defaults
if True:
    sect_name = 'JdFmouth'
    lon0_str = ' -124.6'
    lat0_str = '48.35'
    lon1_str = ' -124.6'
    lat1_str = '48.6'
else:
    sect_name = 'SoGnorth'
    lon0_str = ' -125.4'
    lat0_str = '50.0'
    lon1_str = ' -124.6'
    lat1_str = '50.0'

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v1')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.09.01')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.09.03')
# section specific arguments
parser.add_argument('-sn', '--sect_name', nargs='?', type=str, default=sect_name)
parser.add_argument('-lon0', '--lon0_str', nargs='?', type=str, default=lon0_str)
parser.add_argument('-lat0', '--lat0_str', nargs='?', type=str, default=lat0_str)
parser.add_argument('-lon1', '--lon1_str', nargs='?', type=str, default=lon1_str)
parser.add_argument('-lat1', '--lat1_str', nargs='?', type=str, default=lat1_str)
args = parser.parse_args()

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# make sure the output directory exists
outdir = Ldir['LOo'] + 'extract/'
Lfun.make_dir(outdir)
# name output file
out_fn = (outdir + 'tef_' + Ldir['gtagex'] + '_' +
    args.sect_name + '_' + args.date_string0 + '_' + args.date_string1 + '.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass

dt0 = datetime.strptime(args.date_string0, '%Y.%m.%d')
dt1 = datetime.strptime(args.date_string1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1

# get list of history files to plot
fn_list = Lfun.get_fn_list('hourly', Ldir, args.date_string0, args.date_string1)

# end points
x0 = float(args.lon0_str)
y0 = float(args.lat0_str)
x1 = float(args.lon1_str)
y1 = float(args.lat1_str)

# determine the direction of the section
if (x0==x1) and (y0!=y1):
    sdir = 'NS'
elif (x0!=x1) and (y0==y1):
    sdir = 'EW'
else:
    print('Input points do not form a proper section')
    sdir='bad'
    sys.exit()

print('sect_name = ' + sect_name)
print('sdir = ' + sdir)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)

# we assume a plaid grid, as usual
if sdir == 'NS':
    lon = G['lon_u'][0,:].squeeze()
    lat = G['lat_u'][:,0].squeeze()
elif sdir == 'EW':
    lon = G['lon_v'][0,:].squeeze()
    lat = G['lat_v'][:,0].squeeze()

# we get all 4 i's or j's but only 3 are used
i0, i1, fr = zfun.get_interpolant(np.array([x0]), lon, extrap_nan=True)
if np.isnan(fr):
    print('Bad x point')
    sys.exit()
else:
    ii0 = int(i0)
i0, i1, fr = zfun.get_interpolant(np.array([x1]), lon, extrap_nan=True)
if np.isnan(fr):
    print('Bad x point')
    sys.exit()
else:
    ii1 = int(i1)
j0, j1, fr = zfun.get_interpolant(np.array([y0]), lat, extrap_nan=True)
if np.isnan(fr):
    print('Bad y0 point')
    sys.exit()
else:
    jj0 = int(j0)
j0, j1, fr = zfun.get_interpolant(np.array([y1]), lat, extrap_nan=True)
if np.isnan(fr):
    print('Bad y1 point')
    sys.exit()
else:
    jj1 = int(j1)

# get mask and trim indices
# Note: the mask in G is True on water points
if sdir == 'NS':
    mask = G['mask_u'][jj0:jj1+1, ii0]
    # Note: argmax finds the index of the first True in this case
    igood0 = np.argmax(mask)
    igood1 = np.argmax(mask[::-1])
    # keep one mask point on each end, just to be sure we have a closed section
    Mask = mask[igood0-1:-igood1+1]
    # and change the indices to match.  These will be the indices
    # of the start and end points.
    jj0 = jj0 + igood0 - 1
    jj1 = jj1 - igood1 + 1
    # testing
    Mask_alt = G['mask_u'][jj0:jj1+1, ii0]
    if (Mask==Mask_alt).all():
        print('mask trimmed correctly')
        print('jj0=%d, jj1=%d, ii0=%d' % (jj0, jj1, ii0))
    Lat = lat[jj0:jj1+1]
    Lon = lon[ii0] * np.ones_like(Mask)
elif sdir == 'EW':
    mask = G['mask_v'][jj0, ii0:ii1+1]
    igood0 = np.argmax(mask)
    igood1 = np.argmax(mask[::-1])
    Mask = mask[igood0-1:-igood1+1]
    ii0 = ii0 + igood0 - 1
    ii1 = ii1 - igood1 + 1
    Mask_alt = G['mask_v'][jj0, ii0:ii1+1]
    if (Mask==Mask_alt).all():
        print('mask trimmed correctly')
        print('jj0=%d, ii1=%d, ii1=%d' % (jj0, ii0, ii1))
    Lon = lon[ii0:ii1+1]
    Lat = lat[jj0] * np.ones_like(Mask)
    
# generating some lists
vn_list = []
ds = nc.Dataset(fn)
if False:
    # all 3D variables on the s_rho grid
    for vv in ds.variables:
        vdim = ds.variables[vv].dimensions
        if ( ('ocean_time' in vdim) and ('s_rho' in vdim) ):
            vn_list.append(vv)
else:
    # override
    vn_list.append('salt')
# and some dicts of long names and units
long_name_dict = dict()
units_dict = dict()
for vn in vn_list + ['ocean_time']:
    try:
        long_name_dict[vn] = ds.variables[vn].long_name
    except:
        long_name_dict[vn] = ''
    try:
        units_dict[vn] = ds.variables[vn].units
    except:
        units_dict[vn] = ''
ds.close()
# add custom dict fields
long_name_dict['q'] = 'transport'
units_dict['q'] = 'm3 s-1'
long_name_dict['lon'] = 'longitude'
units_dict['lon'] = 'degrees'
long_name_dict['lat'] = 'latitude'
units_dict['lat'] = 'degrees'

# initialize netcdf output file
NT = len(fn_list)
NX = len(Mask)
NZ = S['N']
foo = nc.Dataset(out_fn, 'w')
foo.createDimension('xi_sect', NX)
foo.createDimension('s_rho', NZ)
foo.createDimension('ocean_time', NT)
foo.createDimension('sdir_str', 2)
for vv in ['ocean_time']:
    v_var = foo.createVariable(vv, float, ('ocean_time',))
    v_var.long_name = long_name_dict[vv]
    v_var.units = units_dict[vv]
for vv in vn_list + ['q']:
    v_var = foo.createVariable(vv, float, ('ocean_time', 's_rho', 'xi_sect'))
    v_var.long_name = long_name_dict[vv]
    v_var.units = units_dict[vv]
for vv in ['lon', 'lat']:
    v_var = foo.createVariable(vv, float, ('xi_sect'))
    v_var.long_name = long_name_dict[vv]
    v_var.units = units_dict[vv]
for vv in ['zeta']:
    v_var = foo.createVariable(vv, float, ('ocean_time', 'xi_sect'))
    v_var.long_name = 'Free Surface Height'
    v_var.units = 'm'

# add static variables
foo['lon'][:] = Lon
foo['lat'][:] = Lat

# add global attributes
foo.sdir = sdir
foo.gtagex = Ldir['gtagex']
foo.date_string0 = args.date_string0
foo.date_string1 = args.date_string1

# extract and save time-dependent fields
count = 0
for fn in fn_list:
    if np.mod(count,24)==0:
        print(' working on %d of %d' % (count, NT))
        sys.stdout.flush()
    ds = nc.Dataset(fn)
    
    # get depth and dz
    if sdir=='NS':
        h = ds['h'][jj0:jj1+1,ii0:ii1+1].squeeze()
        zeta = ds['zeta'][0,jj0:jj1+1,ii0:ii1+1].squeeze()
        z = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(z, axis=0)
        DZ = dz.mean(axis=2)
        dy = G['DY'][jj0:jj1+1,ii0:ii1+1].squeeze()
        DY = dy.mean(axis=1)
        zeta = zeta.mean(axis=1)
    elif sdir=='EW':
        h = ds['h'][jj0:jj1+1,ii0:ii1+1].squeeze()
        zeta = ds['zeta'][0,jj0:jj1+1,ii0:ii1+1].squeeze()
        z = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(z, axis=0)
        DZ = dz.mean(axis=1)
        dy = G['DY'][jj0:jj1+1,ii0:ii1+1].squeeze()
        DY = dy.mean(axis=0)
        zeta = zeta.mean(axis=0)
    # and then create the array of cell areas on the section
    DA = DY.reshape((1, NX)) * DZ
    # then velocity and hence transport
    if sdir=='NS':
        vel = ds['u'][0, :, jj0:jj1+1, ii0].squeeze()
    elif sdir=='EW':
        vel = ds['v'][0, :, jj0, ii0:ii1+1].squeeze()
    q = vel * DA
    
    foo['q'][count, :, :] = q
    foo['zeta'][count, :] = zeta
    foo['ocean_time'][count] = ds['ocean_time'][0]
    
    # save the tracer fields averaged onto this section
    for vn in vn_list:
        if sdir=='NS':
            vv = ds[vn][0,:,jj0:jj1+1,ii0:ii0+2].squeeze()
            vvv = vv.mean(axis=2)
        elif sdir=='EW':
            vv = ds[vn][0,:,jj0:jj0+2,ii0:ii1+1].squeeze()
            vvv = vv.mean(axis=1)
        foo[vn][count,:,:] = vvv
    
    ds.close()
    count += 1

    
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


