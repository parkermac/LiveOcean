"""
Finds the volume of water between TEF sections
"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd

import os
import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart('cas6', 'v2')
import zfun
import zrfun

g_fn = Ldir['grid'] + 'grid.nc'
G = zrfun.get_basic_info(g_fn, only_G = True)

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

# sys.path.append(os.path.abspath(Ldir['parent'] + 'ptools/pgrid'))
# import gfun
# import gfun_plotting as gfp
# Gr = gfun.gstart('cas6')

import tef_fun
from importlib import reload
reload(tef_fun)

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

# choose a pair of sections (Series)
for sn in ['mb1', 'mb2']:
    x0, x1, y0, y1, landward = sect_df.loc[sn,:]
    # also get section orientation
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
        a = [y0, y1]; a.sort()
        y0 = a[0]; y1 = a[1]
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
        a = [x0, x1]; a.sort()
        x0 = a[0]; x1 = a[1]
    ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask = tef_fun.get_inds(x0, x1, y0, y1, G)
    

