"""
Make profiles of salt or density vs. depth at each section.

Goal is to explore time evolution of ds/dx at places like Admiralty Inlet.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import scipy.stats as stats
import pandas as pd

import os
import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun

Ldir = Lfun.Lstart()

indir0 = Ldir['LOo'] + 'tef2/'
# choose the tef extraction to process
#item = Lfun.choose_item(indir0)
year = 2017
item = 'cas6_v3_lo8b_'+str(year)+'.01.01_'+str(year)+'.12.31'
indir0 = indir0 + item + '/'
indir = indir0 + 'extractions/'

sect_list_raw = os.listdir(indir)
sect_list_raw.sort()
sect_list = [item for item in sect_list_raw if ('.nc' in item)]

testing = True

if testing:
    sect_list = ['ai1.nc','ai2.nc','ai3.nc']
else:
    pass
    
outdir = indir0 + 'profiles/'
Lfun.make_dir(outdir)

if testing:
    plt.close('all')

# create z edges and bins
zedges = np.linspace(-100, 0, 26) # original was 1001 used 5001 for Willapa
zbins = zedges[:-1] + np.diff(zedges)/2

for fn in sect_list:
    
    if testing:
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        
    print(fn)
    
    df = pd.DataFrame(index=zbins , columns = range(1,13,1))

    # load fields
    ds = nc.Dataset(indir + fn)
    s = ds['salt'][:]
    mask = s[0,:,:].squeeze().mask
    th = ds['temp'][:]
    z0 = ds['z0'][:]
    ot = ds['ocean_time'][:]
    gtagex = ds.gtagex
    ds0 = ds.date_string0
    ds1 = ds.date_string1
    ds.close()
    
    dt = []
    mo = []
    for t in ot:
        this_dt = Lfun.modtime_to_datetime(t)
        this_mo = this_dt.month
        dt.append(this_dt)
        mo.append(this_mo)
    mo[-1] = 12 # fix last point
    mo = np.array(mo)
    
    if mask.size > 1:
        z0m = z0[~mask].data
    else:
        z0m = z0.flatten()
    
    for mm in (range(1,13,1)):
        i0 = np.argwhere(mo==mm)[0][0]
        i1 = np.argwhere(mo==mm)[-1][0]
        si = s[i0:i1+1,:,:].mean(axis=0)
        if mask.size > 1:
            sm = si[~mask].data
        else:
            sm = si.flatten()
        ss = stats.binned_statistic(z0m, sm, statistic='mean', bins=zedges)
        sb = ss.statistic        
        df.loc[:,mm] = sb
        
    if testing:
        df.plot(ax=ax, title=fn)
        
    df.to_pickle(outdir + fn.replace('.nc','.p'))
    
plt.show()
    

