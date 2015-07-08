"""
This is the main program for making the BIO forcing file.
"""
# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("gridname", type=str, help="cascadia1, etc.")
parser.add_argument("tag", type=str, help="base, etc.")
parser.add_argument("frc", type=str, help="atm, ocn, riv, or tide")
parser.add_argument("run_type", type=str, help="forecast or backfill")
parser.add_argument("date_string", type=str, help="e.g. 2014.02.14")
# and this is an optional input parameter
parser.add_argument("-x", "--ex_name", type=str, help="e.g. lo1")
args = parser.parse_args()
# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['LOogf_f'] = (Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/' + args.frc + '/')
    
# screen output
from datetime import datetime
print('MAIN: frc = ' + args.frc + ', run_type = ' + args.run_type
    + ', date_string = ' + args.date_string)
print('MAIN start time = ' + str(datetime.now()))

# ****************** CASE-SPECIFIC CODE *****************
import shutil
fn_list = []
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/ocn/ocean_clm'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/ocn/ocean_bry'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/ocn/ocean_ini'))
fn_list.append((Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/riv/rivers'))
for fn in fn_list:
    fn_orig = (fn + '.nc')
    fn_for_bio = (fn + '_bio.nc')
    # first get rid of any existing copies of the forcing that
    # might have bio variables, becasue the worker needs them to be absent
    try:
        os.remove(fn_for_bio)
    except OSError:
        pass
        # assume the file did not exist
    # then make copies with new names
    shutil.copyfile(fn_orig,fn_for_bio) 

# ******************************************************* 

# run the code to create the forcing files
#Lfun.run_worker(args.date_string, Ldir)
print('MAIN end time = ' + str(datetime.now()))
