"""
Process a TEF extraction.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'

# choose the mooring extraction to process
print('\n%s\n' % '** Choose TEF extraction file to process **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.nc' in m) and ('tef_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
if False:
    my_npt = int(input('-- Input number -- '))
else:
    my_npt = 0 # for testing
tef_file = m_dict[my_npt]
fn = indir + tef_file

ds = nc.Dataset(fn)

q = ds['q'][:]
salt = ds['salt'][:]
ot = ds['ocean_time'][:]


#ds.close()
