"""
Code for particle tracking.  Meant to be fast and simple, but capable
of tracking backward in time.

Designed for ROMS output with plaid lat, lon grids.

This program is mainly a DRIVER where you supply:
    - which run to work on
    - particle initial locations
    - a list of starting times
    - duration in days to track
    - some other flags for options like using turbulent dispersion.
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
import trackfun_1 as tf1
reload(tf1)

#%% ************ USER INPUT **************************************

# some run specifications
exp_name = 'GCtest'
gtagex = 'cas3_v0_lo6m' # e.g. 'cascadia1_base_lobio1'
ic_name = 'test0' # 'jdf' or 'cr' or etc.
dir_tag = 'forward' # 'forward' or 'reverse'
surface = True # Boolean, True for trap to surface
turb = False # Vertical turbulent dispersion
windage = 0 # a small number >= 0
ndiv = 10 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)

# set run time information (multiple experiments)
#
# always start on a day (no hours)
if Ldir['env'] == 'pm_mac':
    dt_first_day = datetime(2013,1,29)
    number_of_start_days = 1
    days_between_starts = 1
    days_to_track = 1

# set particle initial locations, all numpy arrays
#
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length of pcs00 is however many vertical positions you have at
# each lat, lon

if ic_name in ['test0']:
    lonvec = np.linspace(-122.55, -122.51, 10)
    latvec = np.linspace(47.27, 47.31, 10)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    pcs00 = np.array([-.05])

if len(plon00) != len(plat00):
    print('Problem with length of initial lat, lon vectors')
    sys.exit()
    
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

# info for NetCDF output
name_unit_dict = {'lon':('Longitude','degrees'), 'lat':('Latitude','degrees'),
    'cs':('Fractional Z','Dimensionless'), 'ot':('Ocean Time','Seconds since 1/1/1970 UTC'),
    'z':('Z','m'), 'zeta':('Surface Z','m'), 'zbot':('Bottom Z','m'),
    'salt':('Salinity','Dimensionless'), 'temp':('Potential Temperature','Degrees C'),
    'u':('EW Velocity','meters s-1'), 'v':('NS Velocity','meters s-1'),
    'w':('Vertical Velocity','meters s-1'),
    'Uwind':('EW Wind Velocity','meters s-1'), 'Vwind':('NS Velocity','meters s-1'),
    'h':('Bottom Depth','m')}

#%% step through the experiments (one for each start day)
write_grid_file = True
for idt0 in idt_list:
    
    # split the calculation up into one-day chunks    
    for nd in range(Ldir['days_to_track']):
        
        idt = idt0 + timedelta(days=nd)
        
        # make sure out file list starts at the start of the day
        if nd > 0: 
            fn_last = fn_list[-1]
        fn_list = tf1.get_fn_list_1day(idt, Ldir)
        if nd > 0:
            fn_list = [fn_last] + fn_list
        
        if write_grid_file:
            # write a file of grid info
            dsh = nc4.Dataset(fn_list[0])
            dsg = nc4.Dataset(outdir + 'grid.nc', 'w')
            # lists of variables to process
            dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi']
            vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
            # Copy dimensions
            for dname, the_dim in dsh.dimensions.items():
                if dname in dlist:
                    dsg.createDimension(dname, len(the_dim))
            # Copy variables
            for vn in vn_list2:
                varin = dsh[vn]
                vv = dsg.createVariable(vn, varin.dtype, varin.dimensions)
                vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
                vv[:] = dsh[vn][:]
            dsh.close()
            dsg.close()
            write_grid_file = False
            
        T0 = zrfun.get_basic_info(fn_list[0], only_T=True)
        Ldir['date_string0'] = datetime.strftime(T0['tm'],'%Y.%m.%d')
           
        #%% DO THE TRACKING
        if nd == 0: # first day
            
            plon0 = plon00.copy()
            plat0 = plat00.copy()
            pcs0 = pcs00.copy()
            
            # do the tracking
            tt0 = time.time()    
            P, G, S = tf1.get_tracks(fn_list, plon0, plat0, pcs0,
                                dir_tag, surface, turb,
                                ndiv, windage, trim_loc=True)
            # save the results to NetCDF
            NT, NP = P['lon'].shape
            outname = ('release_' + Ldir['date_string0'] + '_' +
                str(Ldir['days_to_track']) + 'days'+ '.nc')
            print(' - ' + outname)
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
        else: # subsequent days
            # get initial condition
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            # do the tracking
            P, G, S = tf1.get_tracks(fn_list, plon0, plat0, pcs0,
                                dir_tag, surface, turb,
                                ndiv, windage)
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
    print(' - Took %0.1f sec for %s day(s)' %
            (time.time() - tt0, str(Ldir['days_to_track'])))
    print(50*'=')
