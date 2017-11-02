"""
Code to test the subprocess call that is giving me trouble
in python 3.

"""

# setup
import os
import sys
import argparse
from datetime import datetime
# This relative path to alpha is meant to work only when intro()
# is called from the forcing directories (not great coding).
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
# Note: the path "alp" will now also work for the calling function
import Lfun

# set defaults
gridname = 'cascadia1'
tag = 'base'    
cwd = os.getcwd()
icwd = cwd.rfind('/')
frc = cwd[icwd+1:]
run_type = 'forecast'  # backfill or forecast
# Example of date_string is 2015.09.19
date_string = datetime.now().strftime(format='%Y.%m.%d')
ex_name = 'lo1'

# get the dict Ldir
Ldir = Lfun.Lstart(gridname, tag)
Ldir['date_string'] = date_string
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

# add the arguments to Ldir, because some code needs them
Ldir['gridname'] = gridname
Ldir['tag'] = tag
Ldir['frc'] = frc
Ldir['run_type'] = run_type
Ldir['date_string'] = date_string
Ldir['ex_name'] = ex_name

# Make the directory tree for this forcing, if needed. This is redundant
# with what the driver does (except that it clobbers nothing), and is
# only included here so that we can test the python code without using
# the driver.
Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
Ldir['LOogf'] = (Ldir['LOog'] + 'f' + date_string + '/')
Ldir['LOogf_f'] = (Ldir['LOogf'] + frc + '/')
Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')
Lfun.make_dir(Ldir['LOogf'])
Lfun.make_dir(Ldir['LOogf_f'], clean=True)
Lfun.make_dir(Ldir['LOogf_fi'], clean=True)
Lfun.make_dir(Ldir['LOogf_fd'], clean=True)

# pass arguments to a matlab program
import subprocess
func = ("make_forcing_worker(\'" +
    Ldir['gridname'] + "\',\'" +
    Ldir['tag'] + "\',\'" +
    Ldir['date_string'] + "\',\'" +
    Ldir['run_type'] + "\',\'" +
    Ldir['LOogf_f'] + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
print(out.decode())