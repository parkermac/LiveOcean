"""
This calculates hypoxic volume.

"""

# imports

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun
import zfun

import numpy as np
import netCDF4 as nc
import argparse
from datetime import datetime, timedelta
start_time = datetime.now()

# Command line arguments

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2019.07.04')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2019.07.05')
parser.add_argument('-lt', '--list_type', nargs='?', type=str, default='daily')
parser.add_argument('-test', '--testing', default=False, type=boolean_string)
args = parser.parse_args()

# save some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
Ldir['list_type'] = args.list_type
Ldir['date_string0'] = args.date_string0
Ldir['date_string1'] = args.date_string1

# make sure the output directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = outdir00 + 'layer/'
Lfun.make_dir(outdir0)
outdir = (outdir0 + Ldir['gtagex'] + '_' + Ldir['date_string0']
        + '_' + Ldir['date_string1'] + '/')
Lfun.make_dir(outdir, clean=False)

# get list of history files to plot
fn_list = Lfun.get_fn_list(args.list_type, Ldir, args.date_string0, args.date_string1)

# name output file
out_fn = (outdir + 'hypoxic_volume_' + Ldir['list_type'] + '.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass
    
# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list_2d = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list_2d_custom = ['DA']
vn_list_3d_t_custom = ['hyp_dz']

# make some things
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
DA = G['DX'] * G['DY']
ny, nx = DA.shape
h = G['h']
S = zrfun.get_basic_info(fn, only_S=True)
zr, zw = zrfun.get_z(h, 0*h, S)
dzr = np.diff(zw, axis=0)

ds1 = nc.Dataset(fn)
ds2 = nc.Dataset(out_fn, 'w')

# Create dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
        
# Create variables and their attributes
# - first time
vn = 'ocean_time'
varin = ds1[vn]
vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
vv.long_name = varin.long_name
vv.units = varin.units
# - then static 2d fields
for vn in vn_list_2d:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = ds1[vn][:]
# - then static custom fields
for vn in vn_list_2d_custom:
    if vn == 'DA':
        vv = ds2.createVariable(vn, float, ('eta_rho', 'xi_rho'))
        vv.long_name = 'Cell horizontal area'
        vv.units = 'm2'
        vv[:] = DA
#
# - then time-dependent custom 3d fields (processed into 2d)
for vn in vn_list_3d_t_custom:
    if vn == 'hyp_dz':
        vv = ds2.createVariable(vn, float, ('ocean_time', 'eta_rho', 'xi_rho'))
        vv.long_name = 'Thickess of hypoxic water'
        vv.units = 'm'
        vv.time='ocean_time'
# Copy time dependent data
omat = np.nan * np.ones(h.shape)
omat = np.ma.masked_where(G['mask_rho']==0, omat)

tt = 0
NF = len(fn_list)
for fn in fn_list:
    if np.mod(tt,24)==0:
        print(' working on %d of %d' % (tt, NF))
        sys.stdout.flush()
    ds = nc.Dataset(fn)
    
    ds2['ocean_time'][tt] = ds['ocean_time'][0]
        
    for vn in vn_list_3d_t_custom:
        if vn == 'hyp_dz':
            oxy = ds['oxygen'][0,:,:,:]
            dzrm = np.ma.masked_where(oxy > 60, dzr)
            hyp_dz = dzrm.sum(axis=0)
            ds2[vn][tt,:,:] = hyp_dz
        
    tt += 1
    ds.close()

if args.testing:
    sys.path.append(os.path.abspath('../plotting'))
    import pfun
    import matplotlib.pyplot as plt
    plt.close('all')
    fs=16
    plt.rc('font', size=fs)
    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(111)
    cs = ax.pcolormesh(G['lon_psi'],G['lat_psi'],
        hyp_dz[1:-1,1:-1]/h[1:-1,1:-1], vmin=0, vmax=1)
    fig.colorbar(cs)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.contour(G['lon_rho'],G['lat_rho'],h, [100, 200, 2000],
        colors=['darkorange','plum','darkorchid'], linewidths=2, linestyles='solid')
    ax.text(.95,.16,'100 m',color='darkorange',weight='bold',transform=ax.transAxes,ha='right')
    ax.text(.95,.13,'200 m',color='plum',weight='bold',transform=ax.transAxes,ha='right')
    ax.text(.95,.1,'2000 m',color='darkorchid',weight='bold',transform=ax.transAxes,ha='right')
    ax.set_title('Hypoxic Depth Fraction')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-130, -128, -126, -124, -122])
    plt.show()
    plt.rcdefaults()
    
    print(' Saved Variables '.center(50,'='))
    for vn in ds2.variables:
        print(vn)

ds1.close()
ds2.close()

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
    


