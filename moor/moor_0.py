"""
Extract a mooring-like record.

For 24 hour days on my mac this is fast, like 1-2 sec per day, but on
fjord it takes more like 22 sec per day, or 2 hours per year.

This can be run from the linux command line:

This runs with default values (RISE north, three days in September 2015)
python roms_moor0.py

This changes the days to run
python moor_0.py -d0 2015.07.10 -d1 2015.08.01

This gets an IRIS mooring record (on fjord):
python moor_0.py -sn J26A -lon -125.4664 -lat 44.6547 -d0 2013.01.02 -d1 2015.11.01

python moor_0.py -sn J26C -lon -125.4653 -lat 44.6534 -d0 2013.01.02 -d1 2015.11.01

And you can also change the run, the station name and location, etc.
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
import zfun
import netCDF4 as nc

# set defaults
gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio1'
date_string0 = datetime(2015,9,18).strftime(format='%Y.%m.%d')
date_string1 = datetime(2015,9,20).strftime(format='%Y.%m.%d')
list_type = 'backfill' # backfill, low_pass
sta_name = 'RN'
lon_str = '-124.5'
lat_str = '47'

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# optional input arguments
parser.add_argument('-g', '--gridname', nargs='?', const=gridname, type=str, default=gridname)
parser.add_argument('-t', '--tag', nargs='?', const=tag, type=str, default=tag)
parser.add_argument('-x', '--ex_name', nargs='?', const=ex_name, type=str, default=ex_name)
parser.add_argument('-l', '--list_type', nargs='?', const=list_type, type=str, default=list_type)
parser.add_argument('-d0', '--date_string0', nargs='?', const=date_string0, type=str, default=date_string0)
parser.add_argument('-d1', '--date_string1', nargs='?', const=date_string1, type=str, default=date_string1)
parser.add_argument('-sn', '--sta_name', nargs='?', const=sta_name, type=str, default=sta_name)
parser.add_argument('-lon', '--lon_str', nargs='?', const=lon_str, type=str, default=lon_str)
parser.add_argument('-lat', '--lat_str', nargs='?', const=lat_str, type=str, default=lat_str)
args = parser.parse_args()

# get the dict Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
Ldir['list_type'] = args.list_type
Ldir['date_string0'] = args.date_string0
Ldir['date_string1'] = args.date_string1
Ldir['sta_name'] = args.sta_name
Ldir['lon_str'] = args.lon_str
Ldir['lat_str'] = args.lat_str

# make sure the output directory exists
outdir = Ldir['LOo'] + 'moor/'
Lfun.make_dir(outdir)

#%% function definitions
def make_date_list(dt0,dt1,Ldir): # a helpful function
    del_dt = timedelta(1)
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + del_dt
    return date_list

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

#%% set up for the extraction

dt0 = datetime.strptime(Ldir['date_string0'],'%Y.%m.%d') # first day
dt1 = datetime.strptime(Ldir['date_string1'],'%Y.%m.%d') # last day
date_list = make_date_list(dt0,dt1,Ldir)

# target position (-124.5, 47 = RN)
Lon = np.array(float(Ldir['lon_str']))
Lat = np.array(float(Ldir['lat_str']))

# get grid info
indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + date_list[0] + '/'
fn = indir + 'ocean_his_0002.nc'
[G, S] = zfun.get_basic_info(fn, getS=True, getT=False)

# get interpolants for this point
Xi0 = dict(); Yi0 = dict()
Xi1 = dict(); Yi1 = dict()
Aix = dict(); Aiy = dict()
for grd in ['rho', 'u', 'v']:
    xx = G['lon_' + grd][1,:]
    yy = G['lat_' + grd][:,1]
    xi0, xi1, xfr = zfun.get_interpolant(Lon, xx)
    yi0, yi1, yfr = zfun.get_interpolant(Lat, yy)
    Xi0[grd] = xi0
    Yi0[grd] = yi0
    Xi1[grd] = xi1
    Yi1[grd] = yi1
    # create little arrays that are used in the actual interpolation
    Aix[grd] = np.array([1-xfr, xfr]).reshape((1,1,2))
    Aiy[grd] = np.array([1-yfr, yfr]).reshape((1,2))

# generating some lists
v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v']
v1_list = ['ocean_time']
v2_list = []
v3_list_rho = []
v3_list_w = []
ds = nc.Dataset(fn)
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

V = dict()
V_long_name = dict()
V_units = dict()
v_all_list = v0_list + v1_list + v2_list + v3_list_rho + v3_list_w
for vv in v_all_list:
    V[vv] = np.array([])
    try:
        V_long_name[vv] = ds.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = ds.variables[vv].units
    except:
        V_units[vv] = ''

# get static variables
for vv in v0_list:
    xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    vvtemp = ds.variables[vv][yi01, xi01].squeeze()
    V[vv] =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)

ds.close()

#%% extract time-dependent fields

if Ldir['list_type'] == 'backfill':
    # gets 24 hours at a time using MFDataset
    count = 0
    for dd in date_list:
        print('Working on date_list item: ' + dd)
        sys.stdout.flush()
        indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dd + '/'
        fn_list = []
        for hh in range(2,26):
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
        ds = nc.MFDataset(fn_list)
        for vv in v1_list:
            vtemp = ds.variables[vv][:].squeeze()
            V[vv] = np.append(V[vv], vtemp)
        for vv in v2_list:
            xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, yi01, xi01].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            V[vv] = np.append(V[vv], vtemp)
        for vv in v3_list_rho + v3_list_w:
            xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            if count == 0:
                V[vv] = vtemp.T
            else:
                V[vv] = np.concatenate((V[vv], vtemp.T), axis=1)
        # listing of contents, if desired
        if count == 0 and False:
            zfun.ncd(ds)
        count += 1
elif Ldir['list_type'] == 'low_pass':
    # gets one at a time
    count = 0
    for dd in date_list:
        print('Working on date_list item: ' + dd)
        sys.stdout.flush()
        indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dd + '/'
        fn = indir + 'low_passed.nc'
        ds = nc.Dataset(fn)
        for vv in v1_list:
            vtemp = ds.variables[vv][:].squeeze()
            V[vv] = np.append(V[vv], vtemp)
        for vv in v2_list:
            xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, yi01, xi01].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            V[vv] = np.append(V[vv], vtemp)
        for vv in v3_list_rho:
            xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            if count == 0:
                V[vv] = vtemp.reshape((S['N'],1))
            else:
                V[vv] = np.concatenate((V[vv], vtemp.reshape((S['N'],1))), axis=1)
        for vv in v3_list_w:
            xi01, yi01, aix, aiy = get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yi01, xi01].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            if count == 0:
                V[vv] = vtemp.reshape((S['N']+1,1))
            else:
                V[vv] = np.concatenate((V[vv], vtemp.reshape((S['N']+1,1))), axis=1)
        # listing of contents, if desired
        if count == 0 and False:
            zfun.ncd(ds)
        count += 1

ds.close()

# create z_rho and z_w (has to be done after we have V['zeta'])
hh = V['h'][:] * np.ones_like(V['zeta'])
z_rho, z_w = zfun.get_z(hh, V['zeta'][:], S)
V['hh'] = hh
V_long_name['hh'] = 'bottom depth (positive down) as a vector'
V_units['hh'] = 'm'
V['z_rho'] = z_rho
V_long_name['z_rho'] = 'z on rho points (positive up)'
V_units['z_rho'] = 'm'
V['z_w'] = z_w
V_long_name['z_w'] = 'z on w points (positive up)'
V_units['z_w'] = 'm'
v2_list.append('hh')
v3_list_rho.append('z_rho')
v3_list_w.append('z_w')

#%% save the output to NetCDF
out_fn = (outdir +
    Ldir['gtagex'] + '_' +
    Ldir['sta_name'] + '_' +
    Ldir['list_type'] + '_' +
    Ldir['date_string0'] + '_' +
    Ldir['date_string1'] +
    '.nc')
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(out_fn, 'w')

N = S['N']
NT = len(V['ocean_time'][:])

foo.createDimension('scalar', 1)
foo.createDimension('s_rho', N)
foo.createDimension('s_w', N+1)
foo.createDimension('ocean_time', NT)

for vv in v0_list:
    v_var = foo.createVariable(vv, float, ('scalar'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v1_list:
    v_var = foo.createVariable(vv, float, ('ocean_time',))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v2_list:
    v_var = foo.createVariable(vv, float, ('ocean_time',))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v3_list_rho:
    v_var = foo.createVariable(vv, float, ('s_rho', 'ocean_time'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]
for vv in v3_list_w:
    v_var = foo.createVariable(vv, float, ('s_w', 'ocean_time'))
    v_var[:] = V[vv][:]
    v_var.long_name = V_long_name[vv]
    v_var.units = V_units[vv]

foo.close()

