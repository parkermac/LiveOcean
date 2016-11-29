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

# ************** END CASE-SPECIFIC CODE *****************

# run the code to create the forcing files
Lfun.run_worker_post(Ldir)

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))
