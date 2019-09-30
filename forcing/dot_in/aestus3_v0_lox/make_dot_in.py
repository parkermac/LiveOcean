"""
This creates and poulates directories for ROMS runs on gaggle.  It is
designed to work with the "BLANK" version of the .in file,
replacing things like $whatever$ with meaningful values.

"""

import os
import sys
fpth = os.path.abspath('../../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

#import netCDF4 as nc
#import numpy as np
from datetime import datetime, timedelta
fdt = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
fdt_yesterday = fdt - timedelta(1)

print('- dot_in.py creating files for LiveOcean for ' + Ldir['date_string'])

gtag = Ldir['gtag']
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
    dtsec = 60
elif Ldir['blow_ups'] == 1:
    dtsec = 30
elif Ldir['blow_ups'] == 2:
    dtsec = 20
elif Ldir['blow_ups'] == 3:
    dtsec = 10
elif Ldir['blow_ups'] == 4:
    dtsec = 6
elif Ldir['blow_ups'] == 5:
    dtsec = 4
elif Ldir['blow_ups'] == 6:
    dtsec = 2
elif Ldir['blow_ups'] == 7:
    dtsec = 1
else:
    print('Unsupported number of blow ups: %d' % (Ldir['blow_ups']))
    
ndtfast = 20
    
restart_nrrec = '-1' # '-1' for a non-crash restart file, otherwise '1' or '2'
his_interval = 3600 # seconds to define and write to history files
rst_interval = 10 # days between writing to the restart file (e.g. 5)

# which forcings to look for
atm_dir = 'BLANK/' # which atm forcing files to use
ocn_dir = 'ocnA/' # which ocn forcing files to use
riv_dir = 'rivE/' # which riv forcing files to use
tide_dir = 'tideA/' # which tide forcing files to use

#### END USER DEFINED VALUES ####

# DERIVED VALUES

if multi_core:
    if Ldir['np_num'] == 64: # for new mox nodes 2*32=64 2019_02
        ntilei = '8' # number of tiles in I-direction
        ntilej = '8' # number of tiles in J-direction
    elif Ldir['np_num'] == 196:
        ntilei = '14' # number of tiles in I-direction
        ntilej = '14' # number of tiles in J-direction
    else:
        print('Unsupported number of processors: %d' % (Ldir['np_num']))
else:
    ntilei = '1'
    ntilej = '1'

# if np.mod(3600,dtsec) != 0:
#     print('** WARNING: dtsec does not fit evenly into 1 hour **')
if dtsec == int(dtsec):
    dt = str(dtsec) + '.0d0' # a string version of dtsec, for the .in file
else:
    dt = str(dtsec) + 'd0' # a string version of dtsec, for the .in file
ninfo = int(his_interval/dtsec) # how often to write info to the log file (# of time steps)
nhis = int(his_interval/dtsec) # how often to write to the history files
ndefhis = int(nhis) # how often to create new history files
nrst = int(rst_interval*86400/dtsec)
ntimes = int(days_to_run*86400/dtsec)

# file location stuff
date_string = Ldir['date_string']
date_string_yesterday = fdt_yesterday.strftime('%Y.%m.%d')
dstart = str(int(Lfun.datetime_to_modtime(fdt) / 86400.))
f_string = 'f' + date_string
f_string_yesterday = 'f'+ date_string_yesterday
# where forcing files live (fjord, as seen from gaggle)
# NOTE: eventually this should not be hard-wired.
lo_dir = Ldir['parent'] + 'LiveOcean/'
loo_dir = Ldir['parent'] + 'LiveOcean_output/'
grid_dir = Ldir['parent'] + 'LiveOcean_data/grids/' + Ldir['gridname'] + '/'
force_dir = loo_dir + gtag + '/' + f_string + '/'
roms_dir = Ldir['parent'] + 'LiveOcean_roms/'

# hardwired because we don't have netCDF4
nrows = 385 - 2
ncols = 142 - 2

# determine number of layers
s_dict = Lfun.csv_to_dict(grid_dir + 'S_COORDINATE_INFO.csv')
nlayers = str(s_dict['N'])

# the .in file
dot_in_name = 'liveocean.in' # name of the .in file
dot_in_dir00 = Ldir['roms'] + 'output/'
Lfun.make_dir(dot_in_dir00) # make sure it exists
dot_in_dir0 = Ldir['roms'] + 'output/' + gtagex + '/'
Lfun.make_dir(dot_in_dir0) # make sure it exists
dot_in_dir = dot_in_dir0 + f_string +'/'
Lfun.make_dir(dot_in_dir, clean=True) # make sure it exists and is empty

# where to put the output files according to the .in file
out_dir0 = roms_dir + 'output/' + gtagex + '/'
out_dir = out_dir0 + f_string + '/'

if Ldir['start_type'] == 'continuation':
    nrrec = '0' # '-1' for a hot restart?
    ininame = 'ocean_rst.nc' # for a hot perfect restart
    ini_fullname = out_dir0 + f_string_yesterday + '/' + ininame
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini' + '.nc' # could be an ini or history file
    ini_fullname = force_dir + ocn_dir + ininame

# END DERIVED VALUES

## create .in ##########################

f = open('BLANK.in','r')
f2 = open(dot_in_dir + dot_in_name,'w')
in_varlist = ['base_dir','ntilei','ntilej','ntimes','dt','nrrec','ninfo',
    'nhis','dstart','ndefhis','nrst','force_dir','grid_dir','roms_dir',
    'atm_dir','ocn_dir','riv_dir','tide_dir','dot_in_dir',
    'ini_fullname','out_dir','EX_NAME','roms_name',
    'nrows','ncols', 'nlayers', 'ndtfast']
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


