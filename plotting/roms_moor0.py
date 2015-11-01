"""
Extract a mooring-like record.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path: sys.path.append(alp)
import Lfun; reload(Lfun)
import numpy as np
from datetime import datetime, timedelta
import zfun; reload(zfun) # plotting functions
import matfun; reload(matfun) # functions for working with mat files

# set defaults
gridname = 'cascadia1'
tag = 'base'
ex_name = 'lo1'   
date_string0 = datetime(2015,9,18).strftime(format='%Y.%m.%d')
date_string1 = datetime(2015,9,20).strftime(format='%Y.%m.%d')
list_type = 'backfill' # backfill, forecast, low_pass
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
    
def make_date_list(dt0,dt1,Ldir): # a helpful function
    del_dt = timedelta(1)
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + del_dt
    return date_list

dt0 = datetime.strptime(Ldir['date_string0'],'%Y.%m.%d') # first day
dt1 = datetime.strptime(Ldir['date_string1'],'%Y.%m.%d') # last day
date_list = make_date_list(dt0,dt1,Ldir)
    
# target position (-124.5, 47 = RN)
Lon = np.array(float(Ldir['lon_str']))
Lat = np.array(float(Ldir['lat_str']))

# get grid info
indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + date_list[0] + '/'
fn = indir + 'ocean_his_0002.nc'
[G] = zfun.get_basic_info(fn, getS=False, getT=False)

# get interpolants for this point
Xit = dict(); Yit = dict()
Aix = dict(); Aiy = dict()
for grd in ['rho', 'u', 'v']:
    xx = G['lon_' + grd][1,:]
    yy = G['lat_' + grd][:,1]
    xit_list = zfun.get_interpolant(Lon, xx)
    yit_list = zfun.get_interpolant(Lat, yy)
    # this just pulls the interpolant tuple out of the (one-element) list
    # that get_interpolant returns
    Xit[grd] = xit_list[0]; Yit[grd] = yit_list[0]
    # create little arrays that are used in the actual interpolation
    Aix[grd] = np.array([1-Xit[grd][2], Xit[grd][2]]).reshape((1,1,2))           
    Aiy[grd] = np.array([1-Yit[grd][2], Yit[grd][2]]).reshape((1,2))

v1_list = ['ocean_time']
v2_list = ['zeta','sustr','svstr','swrad','lwrad',
    'shflux','latent','sensible']
v3_list = ['temp','salt','u','v']

V = dict()
for vv in v1_list:
    V[vv] = np.array([])
for vv in v2_list:
    V[vv] = np.array([])
for vv in v3_list:
    V[vv] = np.array([])
    
def get_its(ds, vv, Xit, Yit, Aix, Aiy):
    dims = ds.variables[vv].dimensions
    if 'eta_rho' in dims:
        grd = 'rho'
    elif 'eta_u' in dims:
        grd = 'u'
    elif 'eta_v' in dims:
        grd = 'v'
    else:
        print 'grid error!'
    xit = Xit[grd]; yit = Yit[grd]
    aix = Aix[grd]; aiy = Aiy[grd]
    return xit, yit, aix, aiy
    
if Ldir['list_type'] == 'backfill':
    count = 0    
    for dd in date_list:
        print 'Working on date_list item: ' + dd
        sys.stdout.flush()
        indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dd + '/'
        fn_list = []
        for hh in range(2,26):
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
            import netCDF4 as nc
        ds = nc.MFDataset(fn_list)
        for vv in v1_list:
            vtemp = ds.variables[vv][:].squeeze()
            V[vv] = np.append(V[vv], vtemp)
        for vv in v2_list:
            xit, yit, aix, aiy = get_its(ds, vv, Xit, Yit, Aix, Aiy)
            vvtemp = ds.variables[vv][:, yit[:2], xit[:2]].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            V[vv] = np.append(V[vv], vtemp)        
        for vv in v3_list:
            xit, yit, aix, aiy = get_its(ds, vv, Xit, Yit, Aix, Aiy)
            vvtemp = ds.variables[vv][:, :, yit[:2], xit[:2]].squeeze()
            vtemp =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            if count == 0:
                V[vv] = vtemp.T
            else:
                V[vv] = np.concatenate((V[vv], vtemp.T), axis=1)
        ds.close()    
        # listing of contents, if desired
        if count == 0 and False:   
            zfun.ncd(ds)
        count += 1
                   
# save the results
import cPickle as pickle
outname = (outdir +
    Ldir['gtagex'] + '_' +
    Ldir['sta_name'] + '_' +
    Ldir['list_type'] + '_' +
    Ldir['date_string0'] + '_' +
    Ldir['date_string1'] +
    '.p')
pickle.dump( (V, v1_list, v2_list, v3_list, G) , open( outname, 'wb' ) )

