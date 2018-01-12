"""
Code for particle tracking, designed for ROMS output with
plaid lat, lon grids.

This program is a driver where you specify:
- an experiment (ROMS run + release locations + other choices)
- a release or set of releases within that experiment (start day, etc.)

It can be run on its own, or with command line arguments to facilitate
large, automated jobs, for example from the terminal:
python tracker_1.py -dtt 2 -ds 2013.01.30
or in ipython:
run tracker_1.py -dtt 2 -ds 2013.01.30
"""

#%% setup
import numpy as np
from datetime import datetime, timedelta
import time
import netCDF4 as nc4
import argparse
import collections

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
import trackfun_1 as tf1
reload(tf1)
import trackfun_nc as tfnc
reload(tfnc)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-gx', '--gtagex', nargs='?', type=str)
parser.add_argument('-ic', '--ic_name', nargs='?', type=str)
parser.add_argument('-dir', '--dir_tag', nargs='?', type=str)
parser.add_argument('-sfc', '--surface', nargs='?', type=bool)
parser.add_argument('-trb', '--turb', nargs='?', type=bool)
parser.add_argument('-wnd', '--windage', nargs='?', type=float)
parser.add_argument('-ds', '--ds_first_day', nargs='?', type=str)
parser.add_argument('-nsd', '--number_of_start_days', nargs='?', type=int)
parser.add_argument('-dbs', '--days_between_starts', nargs='?', type=int)
parser.add_argument('-dtt', '--days_to_track', nargs='?', type=int)
args = parser.parse_args()
a_dict = args.__dict__ 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ************ USER INPUT **************************************
# The values set in this section are effectively defaults, and
# they will be overridden by the command line arguments.

# (1) specify the experiment
#
exp_name = 'jdf5'
#
if exp_name == 'jdf5':
    gtagex = 'cascadia1_base_lobio5'
    ic_name = 'jdf0'
elif exp_name == 'jdf6':
    gtagex = 'cas3_v0_lo6m'
    ic_name = 'jdf0'
elif exp_name == 'hc5':
    gtagex = 'cascadia1_base_lobio5'
    ic_name = 'hc0'
elif exp_name == 'hc6':
    gtagex = 'cas3_v0_lo6m'
    ic_name = 'hc0'
#
# other experiment choices
dir_tag = 'forward' # 'forward' or 'reverse'
surface = True # Boolean, True to trap to surface
turb = False # Vertical turbulent dispersion
windage = 0 # a small number 0 <= windage << 1 (e.g. 0.03)
# fraction of windspeed added to advection, only for surface=True
if surface == False:
    windage = 0 # override
# additional choices, less likely to change
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration, but still only report fields hourly
#
# modify the experiment name, based on other choices
if dir_tag == 'reverse':
    exp_name = exp_name + '_reverse'
if surface:
    exp_name = exp_name + '_surf'
if turb:
    exp_name = exp_name + '_turb'
if windage > 0:
    exp_name = exp_name + '_wind'
    
# (2) set release information
# 
# You can make multiple releases using:
# number_of_start_days > 1 & days_between_starts
#
if Ldir['env'] == 'pm_mac':
    ds_first_day = '2013.01.29'
    number_of_start_days = 1
    days_between_starts = 1
    days_to_track = 1
    
# routines to set particle initial locations, all numpy arrays
#
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length of pcs00 is however many vertical positions you have at
# each lat, lon (expressed as fraction of depth -1 < pcs < 1)
#
if ic_name == 'hc0': # Hood Canal
    lonvec = np.linspace(-122.65, -122.45, 30)
    latvec = np.linspace(47.2, 47.35, 30)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
    pcs_vec = np.array([-.05])
elif ic_name == 'jdf0': # Mid-Juan de Fuca
    lonvec = np.linspace(-123.85, -123.6, 20)
    latvec = np.linspace(48.2, 48.4, 20)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
    pcs_vec = np.array([-.05])
    
# ********* END USER INPUT *************************************

if args.gtagex != None:
    Ldir['gtagex'] = args.gtagex
else:
    Ldir['gtagex'] = gtagex

# save some things in tr_dict
tr_dict = collections.OrderedDict()
tr_dict['exp_name'] = exp_name
tr_dict['gtagex'] = gtagex
tr_dict['ic_name'] = ic_name
tr_dict['dir_tag'] = dir_tag
tr_dict['surface'] = surface
tr_dict['turb'] = turb
tr_dict['windage'] = windage
tr_dict['ndiv'] = ndiv
tr_dict['ds_first_day'] = ds_first_day
tr_dict['number_of_start_days'] = number_of_start_days
tr_dict['days_between_starts'] = days_between_starts
tr_dict['days_to_track'] = days_to_track

# override using command line arguments
for k in a_dict.keys():
    if a_dict[k] != None:
        tr_dict[k] = a_dict[k]

# make the full IC vectors, which will have equal length
# (one value for each particle)
if len(plon_vec) != len(plat_vec):
    print('Problem with length of initial lat, lon vectors')
    sys.exit()
NSP = len(pcs_vec)
NXYP = len(plon_vec)
plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
plon00 = plon_arr.flatten()
plat00 = plat_arr.flatten()
pcs00 = pcs_arr.flatten()

# make the list of start days (datetimes)
idt_list = []
dt = datetime.strptime(tr_dict['ds_first_day'], '%Y.%m.%d')
for nic in range(number_of_start_days):
    idt_list.append(dt)
    dt = dt + timedelta(days_between_starts)

# make sure the output parent directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

# make the output directory (empty)
outdir1 = tr_dict['exp_name'] + '/'
outdir = outdir0 + outdir1
Lfun.make_dir(outdir, clean=True)
print(50*'*' + '\nWriting to ' + outdir)

# write a csv file of experiment information
Lfun.dict_to_csv(tr_dict, outdir + 'exp_info.csv')

# step through the releases, one for each start day
for idt0 in idt_list:
    
    tt0 = time.time() # monitor integration time
    
    # name the release file by start day
    idt0_str = datetime.strftime(idt0,'%Y.%m.%d')
    outname = ('release_' + idt0_str + '.nc')
    print(' - ' + outname)
    out_fn = outdir + outname
    
    # we do the calculation in one-day segments
    for nd in range(tr_dict['days_to_track']):
        
        # get or replace the history file list for this day
        idt = idt0 + timedelta(days=nd)
        fn_list = tf1.get_fn_list(idt, Ldir)
        # if this is not the first day in the release, we use
        # fn_list_prev to get the first file (hour 0) for this day
        if nd > 0:
            fn_list = [fn_list_prev[-1]] + fn_list
        
        # write the grid file (once per experiment) for plotting
        if idt0 == idt_list[0]:
            g_infile = fn_list[0]
            g_outfile = outdir + 'grid.nc'
            tfnc.write_grid(g_infile, g_outfile)
           
        # DO THE TRACKING
        if nd == 0: # first day
            # set IC
            plon0 = plon00.copy()
            plat0 = plat00.copy()
            pcs0 = pcs00.copy()
            # do the tracking
            P = tf1.get_tracks(fn_list, plon0, plat0, pcs0, tr_dict, trim_loc=True)
            # save the results to NetCDF
            tfnc.start_outfile(out_fn, P)
        else: # subsequent days
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            P = tf1.get_tracks(fn_list, plon0, plat0, pcs0, tr_dict)
            tfnc.append_to_outfile(out_fn, P)
        fn_list_prev = fn_list
            
    print(' - Took %0.1f sec for %s day(s)' %
            (time.time() - tt0, str(tr_dict['days_to_track'])))
    print(50*'=')
