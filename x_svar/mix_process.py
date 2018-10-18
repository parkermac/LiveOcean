"""
Plots layer records for svar extractions.
"""

# setup
import netCDF4 as nc
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

indir0 = Ldir['LOo'] + 'layer/'

# choose the layer extraction to plot
if False:
    item = Lfun.choose_item(indir0)
else:
    item = 'cas4_v2_lo6biom_2017.06.01_2017.07.01'
indir = indir0 + item + '/'
if False:
    infile = Lfun.choose_item(indir, tag='.nc')
    fn = indir + infile
else:
    fn = indir + 'svar_hourly.nc'

ds = nc.Dataset(fn)
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
ot = ds['ocean_time'][:]
aa = pfun.get_aa(ds)
mix = ds['mix'][:]
mask= ds['mask_rho'][:] # 1 on water, 0 on land
DA = ds['DA'][:] # cell horizontal area
x = ds['lon_rho'][:]
y = ds['lat_rho'][:]
ds.close()

days = (ot - ot[0])/86400.
mdays = Lfun.modtime_to_mdate_vec(ot)
mdt = mdates.num2date(mdays) # list of datetimes of data

# neap time limits: centered on 2017.06.16
ii0_neap = 24*14
ii1_neap = 24*17
mmix_neap = mix[ii0_neap:ii1_neap,:,:].mean(axis=0)

# spring time limits: centered on 2017.06.23 Perigean Spring
ii0_spring = 24*21
ii1_spring = 24*24
mmix_spring = mix[ii0_spring:ii1_spring,:,:].mean(axis=0)

# pack everything up and save it
outdir = Ldir['LOo'] + 'svar/'
Lfun.make_dir(outdir)
out_fn = outdir + 'spring_neap_2017.p'

out_dict = {'xp':xp, 'yp':yp, 'x':x, 'y':y, 'aa':aa, 'DA':DA,
    'mmix_neap':mmix_neap, 'mmix_spring':mmix_spring}
    
pickle.dump(out_dict, open(out_fn, 'wb'))
