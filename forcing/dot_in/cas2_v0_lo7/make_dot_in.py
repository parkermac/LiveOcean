"""
This creates and poulates directories for ROMS runs on gaggle.  It is
designed to work with the "BLANK" version of the .in file,
replacing things like $whatever$ with meaningful values.

All variables with a "_a" appended are for the nested grid.

"""

import os
import sys
fpth = os.path.abspath('../../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

Ldir['gridname_a'] = Ldir['gridname'] + 'a'
Ldir['grid_a'] = Ldir['data'] + 'grids/' + Ldir['gridname_a'] + '/'
Ldir['gtag_a'] = Ldir['gridname_a'] + '_' + Ldir['tag']

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
fdt = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
fdt_yesterday = fdt - timedelta(1)

print('- dot_in.py creating files for LiveOcean for ' + Ldir['date_string'])

gtag = Ldir['gtag']
gtag_a = Ldir['gtag_a']
gtagex = gtag + '_' + Ldir['ex_name']
EX_NAME = Ldir['ex_name'].upper()

#### USER DEFINED VALUES ####

# which ROMS code to use
roms_name = 'LO_ROMS'

# account for differences when using biology
do_bio = False

multi_core = True # use more than one core

if Ldir['run_type'] == 'backfill':
    days_to_run = 1.0
elif Ldir['run_type'] == 'forecast':
    days_to_run = float(Ldir['forecast_days'])

# time step in seconds (should fit evenly into 3600 sec)
if Ldir['blow_ups'] == 0:
    dtsec = 45 
    dtsec_a = 15 
elif Ldir['blow_ups'] == 1:
    dtsec = 30
    dtsec_a = 10
else:
    print('Unsupported number of blow ups: %d' % (Ldir['blow_ups']))
restart_nrrec = '-1' # '-1' for a non-crash restart file, otherwise '1' or '2'
his_interval = 3600 # seconds to define and write to history files
rst_interval = 1 # days between writing to the restart file (e.g. 5)

# which forcings to look for
atm_dir = 'atm/' # which atm forcing files to use
ocn_dir = 'ocn1/' # which ocn forcing files to use
riv_dir = 'riv1/' # which riv forcing files to use
tide_dir = 'tide1/' # which tide forcing files to use

#### END USER DEFINED VALUES ####

# DERIVED VALUES

if multi_core:
    if Ldir['np_num'] == 72:
        ntilei = '6' # number of tiles in I-direction (6)
        ntilej = '12' # number of tiles in J-direction (12)
        ntilei_a = '6' # number of tiles in I-direction (6)
        ntilej_a = '12' # number of tiles in J-direction (12)
    elif Ldir['np_num'] == 144:
        ntilei = '8' # number of tiles in I-direction (6)
        ntilej = '18' # number of tiles in J-direction (12)
        ntilei_a = '8' # number of tiles in I-direction (6)
        ntilej_a = '18' # number of tiles in J-direction (12)
    else:
        print('Unsupported number of processors: %d' % (Ldir['np_num']))
else:
    ntilei = '1'
    ntilej = '1'
    ntilei_a = '1'
    ntilej_a = '1'

if np.mod(3600,dtsec) != 0:
    print('** WARNING: dtsec does not fit evenly into 1 hour **')
    
if dtsec == int(dtsec):
    dt = str(dtsec) + '.0d0' # a string version of dtsec, for the .in file
else:
    dt = str(dtsec) + 'd0' # a string version of dtsec, for the .in file

if dtsec_a == int(dtsec_a):
    dt_a = str(dtsec_a) + '.0d0' # a string version of dtsec, for the .in file
else:
    dt_a = str(dtsec_a) + 'd0' # a string version of dtsec, for the .in file

ninfo = int(his_interval/dtsec) # how often to write info to the log file (# of time steps)
nhis = int(his_interval/dtsec) # how often to write to the history files
ndefhis = int(nhis) # how often to create new history files
nrst = int(rst_interval*86400/dtsec)
ntimes = int(days_to_run*86400/dtsec)

ninfo_a = int(his_interval/dtsec_a) # how often to write info to the log file (# of time steps)
nhis_a = int(his_interval/dtsec_a) # how often to write to the history files
ndefhis_a = int(nhis_a) # how often to create new history files
nrst_a = int(rst_interval*86400/dtsec_a)
ntimes_a = int(days_to_run*86400/dtsec_a)

# file location stuff
date_string = Ldir['date_string']
date_string_yesterday = fdt_yesterday.strftime('%Y.%m.%d')
dstart = str(int(Lfun.datetime_to_modtime(fdt) / 86400.))
f_string = 'f' + date_string
f_string_yesterday = 'f'+ date_string_yesterday

# where forcing files live (fjord, as seen from gaggle)
# NOTE: eventually this should not be hard-wired.
if Ldir['env'] == 'pm_mac':
    # testing
    lo_dir = Ldir['LO']
    loo_dir = Ldir['LOo']
    grid_dir = Ldir['data'] + 'grids/' + Ldir['gridname'] + '/'
    grid_dir_a = Ldir['data'] + 'grids/' + Ldir['gridname_a'] + '/'
else:
    lo_dir = '/fjdata1/parker/LiveOcean/'
    loo_dir = '/fjdata1/parker/LiveOcean_output/'
    grid_dir = '/fjdata1/parker/LiveOcean_data/grids/' + Ldir['gridname'] + '/'
    grid_dir_a = '/fjdata1/parker/LiveOcean_data/grids/' + Ldir['gridname_a'] + '/'
force_dir = loo_dir + gtag + '/' + f_string + '/'
force_dir_a = loo_dir + gtag_a + '/' + f_string + '/'
roms_dir = '/pmr1/parker/LiveOcean_roms/'

# determine grid size
gfn = grid_dir + 'grid.nc'
ds = nc.Dataset(gfn)
h = ds['h'][:]
nrows0, ncols0 = h.shape
nrows = nrows0 - 2
ncols = ncols0 - 2
ds.close()
# determine number of layers
s_dict = Lfun.csv_to_dict(grid_dir + 'S_COORDINATE_INFO.csv')
nlayers = str(s_dict['N'])

# determine nested grid size
gfn = grid_dir_a + 'grid.nc'
ds = nc.Dataset(gfn)
h = ds['h'][:]
nrows0, ncols0 = h.shape
nrows_a = nrows0 - 2
ncols_a = ncols0 - 2
ds.close()
# determine number of layers
s_dict = Lfun.csv_to_dict(grid_dir_a + 'S_COORDINATE_INFO.csv')
nlayers_a = str(s_dict['N'])

if do_bio:
    bio_tag = '_bio'
else:
    bio_tag = ''

# the .in file
dot_in_name = 'liveocean.in' # name of the .in file
dot_in_dir0 = Ldir['roms'] + 'output/' + gtagex + '/'
Lfun.make_dir(dot_in_dir0) # make sure it exists
dot_in_dir = dot_in_dir0 + f_string +'/'
Lfun.make_dir(dot_in_dir, clean=True) # make sure it exists and is empty

# where to put the output files according to the .in file
out_dir0 = roms_dir + 'output/' + gtagex + '/'
out_dir = out_dir0 + f_string + '/'

if Ldir['start_type'] == 'continuation':
    nrrec = '0' # '-1' for a hot restart
    #ininame = 'ocean_rst.nc' # for a hot perfect restart
    ininame = 'ocean_his_0025.nc' # for a hot restart
    ininame_a = 'ocean_his_a_0025.nc' # for a hot restart
    ini_fullname = out_dir0 + f_string_yesterday + '/' + ininame
    ini_fullname_a = out_dir0 + f_string_yesterday + '/' + ininame_a
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini' + bio_tag + '.nc' # could be an ini or history file
    ini_fullname = force_dir + ocn_dir + ininame
    ini_fullname_a = force_dir_a + ocn_dir + ininame
nrrec_a = nrrec

# END DERIVED VALUES

## create .in ##########################

f = open('BLANK.in','r')
f2 = open(dot_in_dir + dot_in_name,'w')
in_varlist = ['base_dir','ntilei','ntilej','ntimes','dt','nrrec','ninfo',
    'nhis','dstart','ndefhis','nrst','force_dir','grid_dir','roms_dir',
    'atm_dir','ocn_dir','riv_dir','tide_dir','dot_in_dir',
    'ini_fullname','out_dir','EX_NAME','roms_name','bio_tag',
    'nrows','ncols','nlayers',
    'ntilei_a','ntilej_a','ntimes_a','dt_a', 'nrrec_a', 'ninfo_a',
    'nhis_a','ndefhis_a','nrst_a','force_dir_a','grid_dir_a',
    'ini_fullname_a','nrows_a','ncols_a','nlayers_a']
for line in f:
    for var in in_varlist:
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(eval(var)))
            line = line2
        else:
            line2 = line
    f2.write(line2)
f.close()
f2.close()

## npzd2o_Banas.in ###########

f = open('npzd2o_Banas_BLANK.in','r')
bio_dot_in_name = 'npzd2o_Banas.in'
f3 = open(dot_in_dir + bio_dot_in_name,'w')
in_varlist = ['force_dir','riv_dir','bio_tag']
for line in f:
    for var in in_varlist:
        if '$'+var+'$' in line:
            line2 = line.replace('$'+var+'$', str(eval(var)))
            line = line2
        else:
            line2 = line
    f3.write(line2)
f.close()
f3.close()
