"""
Code for particle tracking.  Meant to be fast and simple, but capable
of tracking backward in time.

Designed for ROMS output with plaid lat, lon grids.

Performance: about 2 minutes for a month on mac.
(forward RK2, dt = 3600).

The design philosophy is that this should be capable of handling both
LiveOcean and older versions, like from PNWTOX, of how files are stored.
It should also work on different platforms (my mac and fjord at this time).

This program is mainly a DRIVER where you supply:
    - which run to work on
    - particle initial locations
    - a list of starting times
    - duration in days to track
    - forward or backward in time
    - some other flags

"""

#%% setup
import numpy as np
from datetime import datetime, timedelta
import time
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
# Note that we override parts of the Ldir logic,
# using Ldir['gtagex'] to identify a run.

from importlib import reload
import zfun
reload(zfun)
import trackfun
reload(trackfun)

#%% ************ USER INPUT **************************************

# some run specifications
gtagex = 'cascadia1_base_lo1' # 'cascadia1_base_lo1' or 'D2005_his'
ic_name = 'deadBirds' # 'jdf' or 'cr'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk2' # 'rk2' or 'rk4'
surface = True # Boolean, True for trap to surface
windage = 0.02 # a small number >= 0
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)

# set run time information (multiple experiments)
# always start on a day (no hours)
if Ldir['parent'] == '/Users/PM5/Documents/':
    # mac version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2015,9,1) 
        number_of_start_days = 1
        days_between_starts = 7
        days_to_track = 7
    elif gtagex == 'D2005_his':
        dt_first_day = datetime(2005,3,17) 
        number_of_start_days = 3
        days_between_starts = 1
        days_to_track = 2
elif Ldir['parent'] == '/data1/parker/':
    # fjord version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2014,11,1) 
        number_of_start_days = 2
        days_between_starts = 3
        days_to_track = 7
        
# set particle initial locations, all numpy arrays
#        
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length ofpcs00 is however many vertical positions you have at
# each lat, lon
if ic_name == 'jdf':
    plon00 = np.array([-124.65])
    plat00 = np.array([48.48])
    pcs00 = np.linspace(-.95,-.05,20)
elif ic_name == 'cr':
    plon00 = np.array([-123.9])
    plat00 = np.array([46.22])
    pcs00 = np.linspace(-.95,-.05,20)
elif ic_name == 'deadBirds':
    lonvec = np.linspace(-127, -123.9, 10)
    latvec = np.linspace(43.5, 49.5, 15)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    pcs00 = np.array([-.05])

if len(plon00) != len(plat00):
    print('Problem with length of initial lat, lon vectors')
    sys.exit()
# ********* END USER INPUT *************************************

#%%

# save some things in Ldir
Ldir['gtagex'] = gtagex
Ldir['ic_name'] = ic_name
Ldir['dir_tag'] = dir_tag
Ldir['method'] = method
Ldir['surface'] = surface
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
plon0 = plon0.flatten()
plat0 = plat0.flatten()
pcs0 = pcs0.flatten()    

# make the list of start days (datetimes)
idt_list = []
dt = dt_first_day 
for nic in range(number_of_start_days):
    idt_list.append(dt)
    dt = dt + timedelta(days_between_starts)

# make sure the output directory exists
outdir = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir)

#%% step through the experiments (one for each start day)
for idt in idt_list:
    
    fn_list = trackfun.get_fn_list(idt, Ldir)
                      
    [T0] = zfun.get_basic_info(fn_list[0], getG=False, getS=False, getT=True)
    [Tend] = zfun.get_basic_info(fn_list[-1], getG=False, getS=False, getT=True)
    Ldir['date_string0'] = datetime.strftime(T0['tm'],'%Y.%m.%d')
    Ldir['date_string1'] = datetime.strftime(Tend['tm'],'%Y.%m.%d')    
                   
    print(50*'*')
    print('Calculating tracks from ' + Ldir['date_string0'] +
          ' to ' + Ldir['date_string1'])
        
    #%% DO THE TRACKING
    tt0 = time.time()
    
    P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0,
                                  dir_tag, method, surface, ndiv, windage) 
                                  
    print(' - Took %0.1f sec for %d days'
          % (time.time() - tt0, Ldir['days_to_track']))
    
    #%% save the results
    outname = (
        Ldir['gtagex'] + '_' +
        Ldir['ic_name'] + '_' +
        Ldir['method'] + '_' +
        'ndiv' + str(Ldir['ndiv']) + '_' +
        Ldir['dir_tag'] + '_' +
        'surface' + str(Ldir['surface']) + '_' +
        'windage' + str(Ldir['windage']) + '_' +
        Ldir['date_string0'] + '_' +
        str(Ldir['days_to_track']) + 'days' +
        '.p')
        
    pickle.dump( (P, G, S, Ldir) , open( outdir + outname, 'wb' ) )
    print('Results saved to:\n' + outname)
    print(50*'*')
