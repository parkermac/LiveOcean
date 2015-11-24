"""
This is the main program for making the RIV forcing file.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime, timedelta
Info = dict()
Info['run_type'] = Ldir['run_type']
if Ldir['run_type'] == 'backfill':   
    dt1 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') + timedelta(3)
    dt0 = datetime.strptime(Ldir['date_string'],'%Y.%m.%d') - timedelta(3)
    Info['datestring_start'] = dt0.strftime('%Y.%m.%d')
    Info['datestring_end'] = dt1.strftime('%Y.%m.%d')
    
# Make a dataframe with info for rivers to get
import Rfun
rdf, rnames = Rfun.get_rivers_dataframe(Ldir)

# Go out and get the data
riv_df_for_matlab, triv_df_for_matlab = Rfun.get_rivers_data(rdf, rnames, Info, Ldir)

# Write the data out for the worker to use
riv_df_for_matlab.to_csv(Ldir['LOogf_fd'] + 'current_river_table.csv')
triv_df_for_matlab.to_csv(Ldir['LOogf_fd'] + 'current_triver_table.csv')

# ******************************************************* 

# run the code to create the forcing files
Lfun.run_worker(Ldir)

print('MAIN end time = ' + str(datetime.now()))
