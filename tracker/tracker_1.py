"""
Code for particle tracking, designed for ROMS output with
plaid lat, lon grids.

This program is a driver where you specify:
- an experiment (ROMS run + release locations + other choices)
- a release within that experiment (start day)
"""

#%% setup
import numpy as np
from datetime import datetime, timedelta
import time
import netCDF4 as nc4

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
import zrfun
import trackfun_1 as tf1
reload(tf1)
import trackfun_nc as tfnc
reload(tfnc)

#%% ************ USER INPUT **************************************

# (1) specify the experiment
#
exp_name = 'hc5'
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
# specify additional experiment info
dir_tag = 'forward' # 'forward' or 'reverse'
surface = False # Boolean, True for trap to surface
turb = False # Vertical turbulent dispersion
windage = 0 # a small number >= 0
#
if dir_tag == 'reverse':
    exp_name = exp_name + '_reverse'
if surface:
    exp_name = exp_name + '_surf'
if turb:
    exp_name = exp_name + '_turb'
if windage > 0:
    exp_name = exp_name + '_wind'
    
# (2) set release information
# - can make multiple releases using:
#   number_of_start_days and days_between_starts
#
# always start on a day (no hours)
if Ldir['env'] == 'pm_mac':
    dt_first_day = datetime(2013,1,29)
    number_of_start_days = 3
    days_between_starts = 1
    days_to_track = 1
#
# set particle initial locations, all numpy arrays
#
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length of pcs00 is however many vertical positions you have at
# each lat, lon
#
if ic_name == 'hc0': # Hood Canal
    lonvec = np.linspace(-122.65, -122.45, 30)
    latvec = np.linspace(47.2, 47.35, 30)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    pcs00 = np.array([-.05])
elif ic_name == 'jdf0': # Mid-Juan de Fuca
    lonvec = np.linspace(-123.85, -123.6, 20)
    latvec = np.linspace(48.2, 48.4, 20)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    pcs00 = np.array([-.05])
#
if len(plon00) != len(plat00):
    print('Problem with length of initial lat, lon vectors')
    sys.exit()
    
# other choices, less likely to change
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration, but still only report fields hourly

# ********* END USER INPUT *************************************

# save some things in Ldir
Ldir['exp_name'] = exp_name
Ldir['gtagex'] = gtagex
Ldir['ic_name'] = ic_name
Ldir['dir_tag'] = dir_tag
Ldir['surface'] = surface
Ldir['turb'] = turb
Ldir['windage'] = windage
Ldir['ndiv'] = ndiv
Ldir['days_to_track'] = days_to_track

# make the full IC vectors, which will have equal length
# (one value for each particle)
NSP = len(pcs00)
NXYP = len(plon00)
plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
plon00 = plon0.flatten()
plat00 = plat0.flatten()
pcs00 = pcs0.flatten()

# make the list of start days (datetimes)
idt_list = []
dt = dt_first_day
for nic in range(number_of_start_days):
    idt_list.append(dt)
    dt = dt + timedelta(days_between_starts)

# make sure the output directory exists
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

# make the output directory
outdir1 = Ldir['exp_name'] + '/'
outdir = outdir0 + outdir1
Lfun.make_dir(outdir, clean=True)
print(50*'*' + '\nWriting to ' + outdir)

# write a csv file of experiment information
import collections
exp_dict = collections.OrderedDict()
exp_list = ['exp_name', 'gtagex', 'ic_name', 'ndiv',
        'dir_tag', 'surface', 'turb', 'windage']
for item in exp_list:
    exp_dict[item] = str(Ldir[item])
Lfun.dict_to_csv(exp_dict, outdir + 'exp_info.csv')

# step through the releases (one for each start day)
for idt0 in idt_list:
    
    idt0_str = datetime.strftime(idt0,'%Y.%m.%d')
    outname = ('release_' + idt0_str + '.nc')
    print(' - ' + outname)
    out_fn = outdir + outname
    
    # we do the calculation in one-day segments
    for nd in range(Ldir['days_to_track']):
        
        # get or replace the history file list for this day
        idt = idt0 + timedelta(days=nd)
        fn_list = tf1.get_fn_list(idt, Ldir)
        # if this is not the first day in the release, we use
        # fn_list_prev to get the first file (hour 0) for this day
        if nd > 0:
            fn_list = fn_list_prev[-1] + fn_list
        
        # write the grid file (once per experiment)
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
            tt0 = time.time()
            P, G, S = tf1.get_tracks(fn_list, plon0, plat0, pcs0,
                                dir_tag, surface, turb,
                                ndiv, windage, trim_loc=True)
            # save the results to NetCDF
            tfnc.start_outfile(out_fn, P)
        else: # subsequent days
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            P, G, S = tf1.get_tracks(fn_list, plon0, plat0, pcs0,
                                dir_tag, surface, turb,
                                ndiv, windage)
            tfnc.append_to_outfile(out_fn, P)
        fn_list_prev = fn_list
            
    print(' - Took %0.1f sec for %s day(s)' %
            (time.time() - tt0, str(Ldir['days_to_track'])))
    print(50*'=')
