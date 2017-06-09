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

from datetime import datetime, timedelta
fdt = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
fdt_yesterday = fdt - timedelta(1)

print('- dot_in.py creating files for LiveOcean for ' + Ldir['date_string'])

#### USER DEFINED VALUES ####

gtag = Ldir['gtag']
gtagex = gtag + '_' + Ldir['ex_name']
EX_NAME = Ldir['ex_name'].upper()

# account for differences when using biology
# NOTE: this is not  robust because it depends on a specific ex_name
do_bio = True

multi_core = True # use more than one core

if Ldir['run_type'] == 'backfill':
    days_to_run = 1.0
elif Ldir['run_type'] == 'forecast':
    days_to_run = float(Ldir['forecast_days'])

dtsec = 30 # time step in seconds INTEGER (should fit evenly into 3600 sec)
restart_nrrec = '-1' # '-1' for a non-crash restart file, otherwise '1' or '2'
his_interval = 3600 # seconds to define and write to history files
rst_interval = 1 # days between writing to the restart file (e.g. 5)

zqt_height = '2.0d0'
zw_height = '10.0d0'

#### END USER DEFINED VALUES ####

# DERIVED VALUES

if multi_core:
    ntilei = '6' # number of tiles in I-direction (6)
    ntilej = '12' # number of tiles in J-direction (12)
else:
    ntilei = '1'
    ntilej = '1'

if float(3600/dtsec) != 3600.0/dtsec:
    print('** WARNING: dtsec does not fit evenly into 1 hour **')
dt = str(dtsec) + '.0d0' # a string version of dtsec, for the .in file
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
lo_dir = '/fjdata1/parker/LiveOcean/'
loo_dir = '/fjdata1/parker/LiveOcean_output/'
grid_dir = '/fjdata1/parker/LiveOcean_data/grids/' + Ldir['gridname'] + '/'
force_dir = loo_dir + gtag + '/' + f_string + '/'
roms_dir = '/pmr1/parker/LiveOcean_roms/'

if do_bio:
    roms_name = 'LO_ROMS'
    bio_tag = '_bio'
else:
    roms_name = 'ROMS'
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

atm_dir = 'atm/' # which atm forcing files to use
ocn_dir = 'ocn/' # which ocn forcing files to use
riv_dir = 'riv/' # which riv forcing files to use
tide_dir = 'tide/' # which tide forcing files to use

if Ldir['start_type'] == 'continuation':
    nrrec = '0' # '-1' for a hot restart
    ininame = 'ocean_rst.nc' # for a hot perfect restart
    #ininame = 'ocean_his_0025.nc' # for a hot restart
    ini_fullname = out_dir0 + f_string_yesterday + '/' + ininame
elif Ldir['start_type'] == 'new':
    nrrec = '0' # '0' for a history or ini file
    ininame = 'ocean_ini' + bio_tag + '.nc' # could be an ini or history file
    ini_fullname = force_dir + ocn_dir + ininame

# END DERIVED VALUES

## create .in ##########################

f = open('BLANK.in','r')
f2 = open(dot_in_dir + dot_in_name,'w')
in_varlist = ['base_dir','ntilei','ntilej','ntimes','dt','nrrec','ninfo',
    'nhis','dstart','ndefhis','nrst','force_dir','grid_dir','roms_dir',
    'atm_dir','ocn_dir','riv_dir','tide_dir','dot_in_dir',
    'zqt_height','zw_height','ini_fullname','out_dir','EX_NAME','roms_name','bio_tag']
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

## npxd2o_Banas.in ###########

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
