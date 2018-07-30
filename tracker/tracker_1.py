"""
Code for particle tracking, designed for ROMS output with
plaid lat, lon grids.

This program is a driver where you specify:
- an experiment (ROMS run + release locations + other choices)
- a release or set of releases within that experiment (start day, etc.)

The main argument you provide is -exp, which is the experiment name, and
is used by experiments.make_ic() to get the gtagex and initial particle
locations.  Other possible commmand line arguments and their defaults
are explained in the argparse section below.

NOTE: To improve usefulness for people other than me, this driver will
first look for user_experiments.py and user_trackfun.py before loading
my versions.  This allows you to create yout own experiments, and modifications
to the tracking (e.g. for diurnal depth behavior) while still being able
to use git pull to update the main code.

It can be run on its own, or with command line arguments to facilitate
large, automated jobs, for example in python:
    
run tracker_1.py -dtt 2 -ds 2013.01.30
run tracker_1.py -3d True -rev True -dtt 5
run tracker_1.py -exp ae1 -3d True -rev True -dtt 5 -nsd 4 -dbs 4 -ds 2013.03.01
run tracker_1.py -exp ae2 -3d True -rev True -dtt 7 -nsd 3 -dbs 4 -ds 2013.03.01

From the terminal or a script you would use "python" instead of "run".

"""

#%% setup
from datetime import datetime, timedelta
import time
import argparse

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
#
if os.path.isfile('user_trackfun.py'):
    import user_trackfun as tf1
else:
    import trackfun_1 as tf1
reload(tf1)
#
import trackfun_nc as tfnc
reload(tfnc)
#
if os.path.isfile('user_experiments.py'):
    import user_experiments as exp
else:
    import experiments as exp
reload(exp)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==
    
# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()

# Set the experiment name
# (details set in experiments.py, or, if it exists, user_experiments.py)
parser.add_argument('-exp', '--exp_name', default='ae0', type=str)

# These are False unless the flags are used with the argument True
# so if you do NOT use these flags the run will be:
# - forward in time
# - trapped to the surface
# - no vertical turbulent diffusion
parser.add_argument('-rev', default=False, type=boolean_string) # reverse time
parser.add_argument('-3d', default=False, type=boolean_string) # do 3d tracking
parser.add_argument('-turb', default=False, type=boolean_string) # include turbulence

# windage = a small number: 0 <= windage << 1 (e.g. 0.03)
# fraction of windspeed added to advection, only for 3d=False
parser.add_argument('-wnd', '--windage', default=0, type=float)

# set the starting day (will be last day for rev=True)
parser.add_argument('-ds', '--ds_first_day', default='2013.03.01', type=str)

# You can make multiple releases using:
# number_of_start_days > 1 & days_between_starts
parser.add_argument('-nsd', '--number_of_start_days', default=1, type=int)
parser.add_argument('-dbs', '--days_between_starts', default=1, type=int)
parser.add_argument('-dtt', '--days_to_track', default=1, type=int)

# number of divisions to make between saves for the integration
# e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
# for the integration, but still only report fields hourly
parser.add_argument('-ndiv', default=1, type=int)

args = parser.parse_args()
TR = args.__dict__ 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# set dependent fields

# overrides      
if TR['3d']:
    TR['windage'] = 0
if TR['rev']:
    TR['turb'] = False

# get experiment info, including initial condition   
TR['gtagex'], ic_name, plon00, plat00, pcs00 = exp.make_ic(TR['exp_name'])

out_name = TR['exp_name']
# modify the output folder name, based on other choices
if TR['rev']:
    out_name += '_reverse'
if TR['3d']:
    out_name += '_3d'
elif not TR['3d']:
    out_name += '_surf'
if TR['turb']:
    out_name += '_turb'
if TR['windage'] > 0:
    out_name += '_wind'

# make the list of start days (datetimes)
idt_list = []
dt = datetime.strptime(TR['ds_first_day'], '%Y.%m.%d')
for nic in range(TR['number_of_start_days']):
    idt_list.append(dt)
    dt = dt + timedelta(TR['days_between_starts'])

# make sure the output parent directory exists
outdir00 = Ldir['LOo']
Lfun.make_dir(outdir00)
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

Ldir['gtagex'] = TR['gtagex']

# make the output directory (empty)
outdir1 = out_name + '/'
outdir = outdir0 + outdir1
Lfun.make_dir(outdir, clean=True)
print(50*'*' + '\nWriting to ' + outdir)

# write a csv file of experiment information
Lfun.dict_to_csv(TR, outdir + 'exp_info.csv')

# calculate total number of hours for NetCDF output
NT_full = (24 * TR['days_to_track']) + 1

# step through the releases, one for each start day
write_grid = True
for idt0 in idt_list:
    
    tt0 = time.time() # monitor integration time
    
    if TR['rev']:
        idt0 = idt0 + timedelta(days=TR['days_to_track']-1)
    else:
        pass
    
    # name the release file by start day
    idt0_str = datetime.strftime(idt0,'%Y.%m.%d')
    outname = ('release_' + idt0_str + '.nc')
    print('-- ' + outname)
    out_fn = outdir + outname
    
    # we do the calculation in one-day segments
    for nd in range(TR['days_to_track']):
        
        # get or replace the history file list for this day
        if TR['rev']:
            idt = idt0 - timedelta(days=nd)
        else:
            idt = idt0 + timedelta(days=nd)
            
        # debugging
        idt_str = datetime.strftime(idt,'%Y.%m.%d')
        print(' - working on ' + idt_str)
            
        fn_list = tf1.get_fn_list(idt, Ldir)
        
        # write the grid file (once per experiment) for plotting
        if write_grid == True:
            g_infile = fn_list[0]
            g_outfile = outdir + 'grid.nc'
            tfnc.write_grid(g_infile, g_outfile)
            write_grid = False
           
        # DO THE TRACKING
        if nd == 0: # first day
            # set IC
            plon0 = plon00.copy()
            plat0 = plat00.copy()
            pcs0 = pcs00.copy()
            # do the tracking
            P = tf1.get_tracks(fn_list, plon0, plat0, pcs0, TR,
                               trim_loc=True)
            # save the results to NetCDF
            if TR['rev']:
                it0 = NT_full - 24*(nd+1) - 1
                it1 = it0 + 25
            else:
                it0 = 0
                it1 = 25

            tfnc.start_outfile(out_fn, P, NT_full, it0, it1)
        else: # subsequent days
            # set IC
            if TR['rev']:
                plon0 = P['lon'][0,:]
                plat0 = P['lat'][0,:]
                pcs0 = P['cs'][0,:]
            else:
                plon0 = P['lon'][-1,:]
                plat0 = P['lat'][-1,:]
                pcs0 = P['cs'][-1,:]
            # do the tracking
            P = tf1.get_tracks(fn_list, plon0, plat0, pcs0, TR)
            
            if TR['rev']:
                it0 = NT_full - 24*(nd+1) - 1
                it1 = it0 + 25
            else:
                it0 = 24*nd
                it1 = it0 + 25

            tfnc.append_to_outfile(out_fn, P, it0, it1)
            
    print(' - Took %0.1f sec for %s day(s)' %
            (time.time() - tt0, str(TR['days_to_track'])))
    print(50*'=')
