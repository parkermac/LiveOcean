"""
This is the main program for making the ATM forcing file.
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# ****************** CASE-SPECIFIC CODE *****************
# NONE
# ************** END CASE-SPECIFIC CODE *****************

# run the code to create the forcing files
Lfun.run_worker(Ldir)

from datetime import datetime
print('MAIN end time = ' + str(datetime.now()))
