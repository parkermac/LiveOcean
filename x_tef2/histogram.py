"""
Make a salinity histogram for a history file.

"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun

Ldir = Lfun.Lstart()
fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f2019.07.04/ocean_his_0001.nc'

ds = nc.Dataset(fn)
salt = ds['salt'][:].squeeze()
ds.close()

G, S, T = zrfun.get_basic_info(fn)
zw = zrfun.get_z(G['h'], 0*G['h'], S, only_w=True)
DZ = np.diff(zw, axis=0)

NR, NC = G['h'].shape

DA = G['DX'] * G['DY']

dv = DZ * DA.reshape((1,NR,NC))

mask = salt.mask

sm = salt[~mask].data
dvm = dv[~mask].data