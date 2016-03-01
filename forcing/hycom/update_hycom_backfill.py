"""
This code runs a sequence of steps to update the HYCOM backfill files.

Steps:
    get
    extrapolate
    combine
    filter
    
"""

# need to set these
gridname = 'cascadia1'
tag = 'base'

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname,tag)

import hfun_backfill as hfb

hfb.get_hycom(Lfun, Ldir, exnum = '91.1')

hfb.extrapolate_hycom(Lfun, Ldir, tag = '91.1')

hfb.combine_hycom(Lfun, Ldir)

hfb.filter_hycom(Lfun, Ldir)




