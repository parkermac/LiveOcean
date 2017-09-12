"""
Code for particle tracking.  Meant to be fast and simple, but capable
of tracking backward in time.

Designed for ROMS output with plaid lat, lon grids.

Recoded the integration scheme around 7/15/2016 to use new version of
get_interpolant.py, and to better handle many particles.

Recoded this driver to split the run up into 1-day chunks, because of
memory issues in large runs.

Recoded 9/12/2017 to use NetCDF for output, instead of pickle.

PERFORMANCE: With the new fast version is took about 12 seconds
for a day of integration for 6-600 particles, and only increased to
19 seconds for 6000 particles.  The old slow version was faster
for < 100 particles, but otherwise became very slow, scaling linearly
with the number of particles.  These tests were all with the cascadia1 grid.

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
import netCDF4 as nc4

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
import zrfun
import trackfun
reload(trackfun)

#%% ************ USER INPUT **************************************

# some run specifications
gtagex = 'cascadia1_base_lobio1' # e.g. 'cascadia1_base_lo1' or 'D2005_his'
ic_name = 'test' # 'jdf' or 'cr' or etc.
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk4' # 'rk2' or 'rk4'
surface = True # Boolean, True for trap to surface
windage = 0.0 # a small number >= 0
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)

# set run time information (multiple experiments)
#
# always start on a day (no hours)
if Ldir['env'] == 'pm_mac':
    # mac version
    if gtagex == 'cascadia1_base_lobio1':
        dt_first_day = datetime(2017,8,5)
        number_of_start_days = 1
        days_between_starts = 1
        days_to_track = 2
elif Ldir['env'] == 'fjord':
    # fjord version
    pass

# set particle initial locations, all numpy arrays
#
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length of pcs00 is however many vertical positions you have at
# each lat, lon

if ic_name in ['test']:
    lonvec = np.linspace(-127, -123.9, 20)
    latvec = np.linspace(43.5, 49.5, 30)
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
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

# info for NetCDF output
name_unit_dict = {'lon':('Longitude','degrees'), 'lat':('Latitude','degrees'),
    'cs':('Fractional Z','Dimensionless'), 'ot':('Ocean Time','Seconds since 1/1/1970'),
    'z':('Z','m'), 'zeta':('Surface Z','m'), 'zbot':('Bottom Z','m'),
    'salt':('Salinity','Dimensionless'), 'temp':('Potential Temperature','Degrees C'),
    'u':('EW Velocity','meters s-1'), 'v':('NS Velocity','meters s-1'),
    'w':('Vertical Velocity','meters s-1'),
    'Uwind':('EW Wind Velocity','meters s-1'), 'Vwind':('NS Velocity','meters s-1'),
    'h':('Bottom Depth','m')}

#%% step through the experiments (one for each start day)
for idt0 in idt_list:

    # split the calculation up into one-day chunks    
    for nd in range(Ldir['days_to_track']):
        
        idt = idt0 + timedelta(days=nd)
        
        # make sure out file list starts at the start of the day
        if nd > 0: 
            fn_last = fn_list[-1]
        fn_list = trackfun.get_fn_list_1day(idt, Ldir)        
        if nd > 0:
            fn_list = [fn_last] + fn_list
            
        T0 = zrfun.get_basic_info(fn_list[0], only_T=True)
        Tend = zrfun.get_basic_info(fn_list[-1], only_T=True)
        Ldir['date_string0'] = datetime.strftime(T0['tm'],'%Y.%m.%d')
        Ldir['date_string1'] = datetime.strftime(Tend['tm'],'%Y.%m.%d')
    
        print(50*'*')
        print('Calculating tracks from ' + Ldir['date_string0'] +
              ' to ' + Ldir['date_string1'])
           
        #%% DO THE TRACKING
        if nd == 0: # first day
            # make the output directory            
            outdir = (outdir0 +
                Ldir['gtagex'] + '_' +
                Ldir['ic_name'] + '_' +
                Ldir['method'] + '_' +
                'ndiv' + str(Ldir['ndiv']) + '_' +
                Ldir['dir_tag'] + '_' +
                'surface' + str(Ldir['surface']) + '_' +
                'windage' + str(Ldir['windage']) + '_' +
                Ldir['date_string0'] + '_' +
                str(Ldir['days_to_track']) + 'days/')
            Lfun.make_dir(outdir, clean=True)
            # do the tracking
            tt0 = time.time()    
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0,
                                          dir_tag, method, surface, ndiv, windage)    
            print(' - Took %0.1f sec for 1 day' % (time.time() - tt0))        
            # save the results to NetCDF
            NT, NP = P['lon'].shape
            outname = 'tracks.nc'
            out_fn = outdir + outname
            ds = nc4.Dataset(out_fn, 'w')
            ds.createDimension('Time', None)
            ds.createDimension('Particle', NP)
            # Copy variables
            for vn in P.keys():
                varin = P[vn]
                if vn == 'ot':
                    vv = ds.createVariable(vn, float, ('Time'))
                else:
                    vv = ds.createVariable(vn, float, ('Time', 'Particle'))
                vv.long_name = name_unit_dict[vn][0]
                vv.units = name_unit_dict[vn][1]
                vv[:] = P[vn][:]
            ds.close()
            print('Results saved to:\n' + outname)
            print(50*'=')
        else: # subsequent days
            tt0 = time.time()
            # get initial condition
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            # do the tracking
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0,
                                          dir_tag, method, surface, ndiv, windage)    
            print(' - Took %0.1f sec for 1 day' % (time.time() - tt0))
            # save the results
            ds = nc4.Dataset(out_fn, 'a')
            NTx, NPx = ds['lon'][:].shape
            for vn in P.keys():
                varin = P[vn]
                if vn == 'ot':
                    ds[vn][NTx:] = P[vn][1:]
                else:
                    ds[vn][NTx:,:] = P[vn][1:,:]
            ds.close()
            print(50*'-')
