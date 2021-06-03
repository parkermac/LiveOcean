"""
Extract fields at a number of sections which may be used later for TEF analysis
of transport and transport-weighted properties.

All input parameters specified at the command line, so this can be run in the background
because it can take a few hours.  Use "-sn all" to get all sections.

Takes about 15 hours for 39 cas6 sections, per year.

"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
start_time = datetime.now()

import zrfun
import tef_fun
from importlib import reload
reload(tef_fun)

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2019.07.04')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2019.07.05')
parser.add_argument('-sn', '--sect_name', nargs='?', type=str, default='ai1')
# this optional argument allows us to access runs kept in non-standard places
# e.g. -rundir /pmr2/darr/LiveOcean_roms/ for cascadia1_base_lobio5
parser.add_argument('-rundir', '--run_directory', nargs='?', type=str, default='')
# section specific arguments
args = parser.parse_args()

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
# overide model run location
if len(args.run_directory) > 0:
    Ldir['roms'] = args.run_directory
else:
    pass
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
# get time limits
ds0 = args.date_string0; ds1 = args.date_string1
Ldir['date_string0'] = ds0; Ldir['date_string1'] = ds1
dt0 = datetime.strptime(ds0, '%Y.%m.%d'); dt1 = datetime.strptime(ds1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1

print('Working on:')
print(Ldir['gtagex'] + '_' + ds0 + '_' + ds1 +'\n')

# make sure the output directories exist
outdir000 = Ldir['LOo']
Lfun.make_dir(outdir000)
outdir00 = outdir000 + 'tef2/'
Lfun.make_dir(outdir00)
outdir0 = (outdir00 + Ldir['gtagex'] + '_' + ds0 + '_' + ds1 + '/')
Lfun.make_dir(outdir0)
outdir = outdir0 + 'extractions/'
Lfun.make_dir(outdir)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()
# initialize a dictionary of info for each section
sect_info = dict()
# select which sections to extract
if args.sect_name == 'all':
    # full list
    sect_list = [item for item in sect_df.index]
else: # single item
    if args.sect_name in sect_df.index:
        sect_list = [args.sect_name]
    else:
        print('That section is not available')
        sys.exit()

# get list of history files to process
if 'LO_roms' in args.run_directory:
    LO_version = True
else:
    LO_version = False
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1, LO_version=LO_version)
NT = len(fn_list)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
NZ = S['N']

print('\nGetting section definitions and indices')
for sect_name in sect_list:
    print(sect_name)
    # name output file
    out_fn = (outdir + sect_name + '.nc')
    # get section lat, lon, and other info
    x0, x1, y0, y1 = sect_df.loc[sect_name,:]
    # get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    NX = len(Mask) # this the length of the section, not specific to x or y
    # save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, NT, NX, NZ, out_fn)
    # initialize a netcdf file for this section, and get list of variables to gather
    vn_list = tef_fun.start_netcdf(fn, out_fn, NT, NX, NZ, Lon, Lat, Ldir)
    # note that this function deletes the existing out_fn, and also creates the list
    # of variables to extract.

# extract and save time-dependent fields
count = 0
print('\nStarting extraction of fields')
print(vn_list)
for fn in fn_list:
    if np.mod(count,24)==0:
        print('  working on %d of %d' % (count, NT))
        sys.stdout.flush()
    ds = nc.Dataset(fn)
    # loop over all sections
    for sect_name in sect_list:
        sinfo = sect_info[sect_name]
        # this is where we add the data from this history file
        # to all of the sections, each defined by sinfo
        tef_fun.add_fields(ds, count, vn_list, G, S, sinfo)
    ds.close()
    count += 1

# finale
import collections
result_dict = collections.OrderedDict()
result_dict['outdir'] = outdir
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


