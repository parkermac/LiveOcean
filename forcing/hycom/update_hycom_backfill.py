"""
This code runs a sequence of steps to update the HYCOM backfill files.

Steps:
    get
    extrapolate
    combine
    filter
    
"""

# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

import hfun_backfill as hfb
hfb.get_hycom(exnum = '91.1')
hfb.extrapolate_hycom(tag = '91.1')
hfb.combine_hycom()
hfb.filter_hycom()



