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

import subprocess

Ldir['indir'] = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + Ldir['date_string'] + '/'

# make a list of all history files in the directory
his_list_raw = os.listdir(Ldir['indir'])
his_list = [hh for hh in his_list_raw if 'ocean_his' in hh]
his_list.sort()
h_list = [int(his[-7:-3]) for his in his_list]

nf = 15

ii = 0
h0 = h_list[ii*nf]
h1 = h_list[(ii+1)*nf - 1]
print('h0 = %d, h1 =%d' % (h0, h1))
Ldir['h0'] = str(h0)
Ldir['h1'] = str(h1)
# run the code to create the forcing files
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
run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
#subprocess.run(run_cmd)

while h1 < h_list[-1]:
    ii += 1
    h0 = h1 + 1
    try:
        h1 = h_list[(ii+1)*nf - 1]
    except IndexError:
        h1 = h_list[-1]
    print('h0 = %d, h1 =%d' % (h0, h1))
    Ldir['h0'] = str(h0)
    Ldir['h1'] = str(h1)
    # run the code to create the forcing files
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
    run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
    #subprocess.run(run_cmd)

# ************** END CASE-SPECIFIC CODE *****************

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))
