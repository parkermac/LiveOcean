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
from importlib import reload
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
# Note that we override parts of the Ldir logic,
# using Ldir['gtagex'] to identify a run.
import numpy as np
from datetime import datetime, timedelta
import zfun
reload(zfun)
import trackfun
reload(trackfun)
import time

# ************ USER INPUT **************************************

gtagex = 'cascadia1_base_lo1' # 'cascadia1_base_lo1' or 'D2005_his'
ic_name = 'test' # 'jdf' or 'cr'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk2' # 'rk2' or 'rk4'
surface = True
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)

if Ldir['parent'] == '/Users/PM5/Documents/':
    # mac version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2015,9,1) # always start on a day (no hours)
        number_of_start_days = 1 #3
        days_between_starts = 7
        days_to_track = 7 #7
    elif gtagex == 'D2005_his':
        dt_first_day = datetime(2005,3,17) # always start on a day (no hours)
        number_of_start_days = 3
        days_between_starts = 1
        days_to_track = 2
elif Ldir['parent'] == '/data1/parker/':
    # fjord version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2014,1,1) # always start on a day (no hours)
        number_of_start_days = 3 #3
        days_between_starts = 7
        days_to_track = 7 #7
        
# set particle initial locations
if ic_name == 'jdf':
    plon00 = np.array([-124.65]) #np.array([-124.51, -124.50, -124.49])
    plat00 = np.array([48.48]) #np.array([48.41, 48.45, 48.49])
    pcs00 = np.linspace(-.95,-.05,20)
    NSP = len(pcs00)
    NXYP = len(plon00)
    plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
elif ic_name == 'cr':
    plon00 = np.array([-123.9])
    plat00 = np.array([46.22])
    pcs00 = np.linspace(-.95,-.05,20)
    NSP = len(pcs00)
    NXYP = len(plon00)
    plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
elif ic_name == 'test':
    lonvec = np.linspace(-125.5, -124.5, 5)
    latvec = np.linspace(46, 47, 5)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    #pcs00 = np.linspace(-.95,-.05,3)
    pcs00 = np.array([-.05])
    NSP = len(pcs00)
    NXYP = len(plon00)
    plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs0 = np.ones((NXYP,NSP)) * pcs00.reshape(1,NSP)
    
# ********* END USER INPUT *************************************

plon0 = plon0.flatten()
plat0 = plat0.flatten()
pcs0 = pcs0.flatten()

idt_list = []
dt = dt_first_day 
for nic in range(number_of_start_days):
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
        for nday in range(days_to_track):
            fdt = idt + timedelta(nday)
            date_list.append(fdt.strftime('%Y.%m.%d'))
        fn_list = []
        for dd in date_list:
            indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                    '/f' + dd + '/')
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
    Ldir['method'] = method
    Ldir['ndiv'] = str(ndiv)
                
    print(50*'*')
    print('Calculating tracks for ' + Ldir['date_string0'] +
          ' to ' + Ldir['date_string1'])
        
    # DO THE TRACKING
    tt0 = time.time()
    P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, delta_t,
                                  dir_tag, method, surface, ndiv)             
    print(' - Took %0.1f sec for %d days'
          % (time.time() - tt0, round((dt1 - dt0).total_seconds()/86400.)))
    
    #%% save the results
    import pickle
    outname = (
        Ldir['gtagex'] + '_' +
        Ldir['ic_name'] + '_' +
        Ldir['method'] + '_' +
        Ldir['dir_tag'] + '_' +
        'ndiv' + Ldir['ndiv'] + '_' +
        Ldir['date_string0'] + '_' +
        str(days_to_track) + 'days' +
        '.p')
    pickle.dump( (P, G, S, Ldir) , open( outdir + outname, 'wb' ) )
    print('Results saved to:\n' + outname)
    print(50*'*')
