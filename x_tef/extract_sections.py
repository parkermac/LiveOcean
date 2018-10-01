"""
Extract fields at a number of sections which may be used later for TEF analysis
of transport and transport-weighted properties.

Performance - for all sections (primary ones, not side channels) with just salt
this takes about 1.5 minutes per day for cas4 on boiler, so about 10 hours for a year.

To do a year of a shelf section (two variables, for Sam) in cascadia1_base_lobio5
takes about an hour and a half on fjord.  So it should take 7.5 hours for 2013-2017,
and about a full day for all three sections.  The fields are about 900 MB per section
per year.

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
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.09.01')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.09.03')
# section specific arguments
args = parser.parse_args()

# Get Ldir
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
Ldir['date_string0'] = args.date_string0
Ldir['date_string1'] = args.date_string1

# make sure the output directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = outdir00 + 'tef/'
Lfun.make_dir(outdir0)
outdir = (outdir0 + Ldir['gtagex'] + '_' + Ldir['date_string0']
        + '_' + Ldir['date_string1'] + '/')
Lfun.make_dir(outdir, clean=False)

dt0 = datetime.strptime(args.date_string0, '%Y.%m.%d')
dt1 = datetime.strptime(args.date_string1, '%Y.%m.%d')
ndays = (dt1-dt0).days + 1

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, args.date_string0, args.date_string1)
NT = len(fn_list)

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
NZ = S['N']

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# initialize a dictionary of info for each section
sect_info = dict()

# select which sections to extract
if False: # limited list
    sect_list = [item for item in sect_df.index if item in ['shelf_45','shelf_46','shelf_47']]
else: # full list
    sect_list = [item for item in sect_df.index]

print('\nGetting section definitions and indices')
for sect_name in sect_list:
    print(sect_name)
    # name output file
    out_fn = (outdir + sect_name + '.nc')
    # get section lat, lon, and other info
    x0, x1, y0, y1, landward = sect_df.loc[sect_name,:]
    # get indices for this section
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    NX = len(Mask) # this the length of the section, not specific to x or y
    # save some things for later use
    sect_info[sect_name] = (ii0, ii1, jj0, jj1, sdir, landward, NT, NX, NZ, out_fn)
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


