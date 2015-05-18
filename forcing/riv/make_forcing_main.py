"""
This is the main program for making the RIV forcing file.
"""
# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("frc", type=str, help="atm, ocn, riv, or tide")
parser.add_argument("run_type", type=str, help="forecast or backfill")
parser.add_argument("date_string", type=str, help="e.g. 2014.02.14")
args = parser.parse_args()

# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
Ldir['LOogf_f'] = (Ldir['LOo'] + Ldir['gtag'] +
    '/f' + args.date_string + '/' + args.frc + '/')
    
# screen output
from datetime import datetime
print('MAIN: frc = ' + args.frc + ', run_type = ' + args.run_type
    + ', date_string = ' + args.date_string)
print('MAIN start time = ' + str(datetime.now()))

# ****************** CASE-SPECIFIC CODE *****************

from datetime import timedelta
Info = dict()
Info['run_type'] = args.run_type
if args.run_type == 'backfill':   
    dt1 = datetime.strptime(args.date_string,'%Y.%m.%d') + timedelta(3)
    dt0 = datetime.strptime(args.date_string,'%Y.%m.%d') - timedelta(3)
    Info['datestring_start'] = dt0.strftime('%Y.%m.%d')
    Info['datestring_end'] = dt1.strftime('%Y.%m.%d')
    
# Make a dataframe with info for rivers to get
import Rfun; reload(Rfun)
rdf, rnames = Rfun.get_rivers_dataframe(Ldir)

# Go out and get the data
riv_df_for_matlab, triv_df_for_matlab = Rfun.get_rivers_data(rdf, rnames, Info, Ldir)

# Write the data out for the worker to use
Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')
riv_df_for_matlab.to_csv(Ldir['LOogf_fd'] + 'current_river_table.csv')
triv_df_for_matlab.to_csv(Ldir['LOogf_fd'] + 'current_triver_table.csv')

# ******************************************************* 

# run the code to create the forcing files
Lfun.run_worker(args.date_string, Ldir)
print('MAIN end time = ' + str(datetime.now()))
