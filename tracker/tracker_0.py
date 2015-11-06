"""
Code for particle tracking.  Meant to be fast and simple, but capable
of tracking backward in time.

Designed for ROMS output with plaid lat, lon grids.

Performance: about 2 minutes for a month on mac.
(forward half-step, dt = 3600).

The design philosophy is that this should be capable of handling both
LiveOcean and older versions, like from PNWTOX, of how files are stored.
It should also work on platform (my mac and fjord at this time).

This program is mainly a DRIVER where you supply:
    - which run to work on
    - particle initial locations
    - a list of starting times
    - duration in days to track
    - forward or backward in time

"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path: sys.path.append(alp)
import Lfun; reload(Lfun)
# get Ldir
Ldir = Lfun.Lstart()
# Note that we override parts of the Ldir logic,
# using Ldir['gtagex'] to identify a run.
import numpy as np
from datetime import datetime, timedelta
import zfun; reload(zfun) # plotting functions
import trackfun; reload(trackfun)
import time

# ************ USER INPUT **************************************

gtagex = 'D2005_his'
ic_name = 'cr'
dir_tag = 'reverse' # forward or reverse

if gtagex == 'D2005_his':
    dt = datetime(2005,3,17) # always start on a day (no hours)
    days_between_starts = 1
    days_to_track = 2
elif gtagex == 'cascadia1_base_lo1':
    dt = datetime(2015,9,1) # always start on a day (no hours)
    days_between_starts = 1
    days_to_track = 2

# set particle initial locations
if ic_name == 'jdf':
    NP = 20
    plon0 = -124.7 *np.ones(NP)
    plat0 = 48.48 * np.ones(NP)
    pcs0 = np.linspace(-.95,-.05,NP)
elif ic_name == 'cr':
    NP = 20
    plon0 = -123.9 *np.ones(NP)
    plat0 = 46.22 * np.ones(NP)
    pcs0 = np.linspace(-.95,-.05,NP)
    
# ********* END USER INPUT *************************************

idt_list = []    
for nic in range(3):
    idt_list.append(dt)
    dt = dt + timedelta(days_between_starts)

# make sure the output directory exists
outdir = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir)

for idt in idt_list:
    
    Ldir['gtagex'] = gtagex
    if Ldir['gtagex'] == 'cascadia1_base_lo1':
        # LiveOcean version
        Ldir['gtagex'] = 'cascadia1_base_lo1'
        # make the list of input history files
        date_list = []
        for dt in idt_list:
            date_list.append(dt.strftime('%Y.%m.%d'))
            dt = dt + timedelta(days_to_track)
        fn_list = []
        for dd in date_list:
            indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dd + '/'
            for hh in range(2,26):
                hhhh = ('0000' + str(hh))[-4:]
                fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
                
    elif Ldir['gtagex'] == 'D2005_his':
        # Other ROMS runs version
        indir = '/Users/PM5/Documents/roms/output/' + Ldir['gtagex'] + '/'
        save_num_list = range(1,365*24)
        save_dt_list = []
        dt00 = datetime(2005,1,1)
        save_dt_list.append(dt00)
        for sn in save_num_list:
            save_dt_list.append(dt00 + timedelta(hours=sn))
        # keys of this dict are datetimes, and values are history numbers
        save_dt_num_dict = dict(zip(save_dt_list,save_num_list))
        fn_list = []
        for hh in range(days_to_track*24 + 1):
            hh = save_dt_num_dict[idt + timedelta(hours=hh)]
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
        [T0] = zfun.get_basic_info(fn_list[0], getG=False, getS=False, getT=True)
        [T1] = zfun.get_basic_info(fn_list[-1], getG=False, getS=False, getT=True)
        dt0 = T0['tm']
        dt1 = T1['tm']
        Ldir['date_string0'] = datetime.strftime(dt0,'%Y.%m.%d')
        Ldir['date_string1'] = datetime.strftime(dt1,'%Y.%m.%d')
    
    # time step in seconds
    [T00] = zfun.get_basic_info(fn_list[0], getG=False, getS=False)
    [T11] = zfun.get_basic_info(fn_list[1], getG=False, getS=False)
    delta_t = T11['ocean_time'] - T00['ocean_time']
    
    Ldir['ic_name'] = ic_name
    Ldir['dir_tag'] = dir_tag
                
    print(50*'*')
    print('Calculating tracks for ' + Ldir['date_string0'] + ' to ' + Ldir['date_string1'])
        
    # DO THE TRACKING
    tt0 = time.time()
    Plon, Plat, Pcs, Pot, G = trackfun.get_tracks(fn_list, plon0, plat0, pcs0,
        delta_t, dir_tag)              
    print(' - Took %0.1f sec for %d days' % (time.time() - tt0, (dt1 - dt0).days))
    
    # save the results
    import cPickle as pickle
    outname = (
        Ldir['gtagex'] + '_' +
        Ldir['ic_name'] + '_' +
        Ldir['dir_tag'] + '_' +
        Ldir['date_string0'] + '_' +
        Ldir['date_string1'] +
        '.p')
    pickle.dump( (Plon, Plat, Pcs, G, Ldir) , open( outdir + outname, 'wb' ) )
    print('Results saved to:\n' + outname)
    print(50*'*')
