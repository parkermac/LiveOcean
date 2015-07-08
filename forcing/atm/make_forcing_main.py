"""
This is the main program for making the ATM forcing file.
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
# NONE
# ******************************************************* 

# run the code to create the forcing files
Lfun.run_worker(args.date_string, Ldir)
print('MAIN end time = ' + str(datetime.now()))
