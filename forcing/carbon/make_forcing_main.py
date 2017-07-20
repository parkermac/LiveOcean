"""
This is the main program for making the CARBON variable additions
to the ROMS history files.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
# Most everything is done in the matlab worker, but we do
# use a different version of Lfun.run_worker() because it allows
# us to tell the worker where the history files are.

Ldir['indir'] = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + Ldir['date_string'] + '/'

Ldir['h0'] = str(2)
Ldir['h1'] = str(4)
# run the code to create the forcing files
import subprocess

func = ("make_forcing_worker(\'" +
    Ldir['gridname'] + "\',\'" +
    Ldir['tag'] + "\',\'" +
    Ldir['date_string'] + "\',\'" +
    Ldir['run_type'] + "\',\'" +
    Ldir['indir'] + "\',\'" +
    Ldir['h0'] + "\',\'" +
    Ldir['h1'] + "\',\'" +
    Ldir['LOogf_f'] + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
subprocess.run(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


#Lfun.run_worker_post(Ldir)

Ldir['h0'] = str(5)
Ldir['h1'] = str(7)
func = ("make_forcing_worker(\'" +
    Ldir['gridname'] + "\',\'" +
    Ldir['tag'] + "\',\'" +
    Ldir['date_string'] + "\',\'" +
    Ldir['run_type'] + "\',\'" +
    Ldir['indir'] + "\',\'" +
    Ldir['h0'] + "\',\'" +
    Ldir['h1'] + "\',\'" +
    Ldir['LOogf_f'] + "\')")
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
subprocess.run(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# run the code to create the forcing files
#Lfun.run_worker_post(Ldir)

# ************** END CASE-SPECIFIC CODE *****************

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))
