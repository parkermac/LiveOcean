"""
Plot time series of tracer in all volumes, vs. time, from
any of the flux_engine experiments.

"""

# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'
import zfun

import tef_fun
import flux_fun


# Input directory
indir0 = Ldir['LOo'] + 'tef/'
outdir = indir0 + 'misc_figs_cas6/'
indir = indir0 + 'flux_engine/' + Ldir['gtagex'] + '/'
voldir = indir0 + 'volumes_' + Ldir['gridname'] + '/'

# load a Series of the volumes of each segment, created by flux_get_vol.py
v_df = pd.read_pickle(voldir + 'volumes.p')
V = flux_fun.get_V(v_df)

infile = 'IC_SoG_2019_winter.p'
#infile = 'IC_SouthSound_2019_winter.p'

# infile is like "'IC_Whidbey_2019_winter.p'"
source = infile.split('_')[0] + '_' + infile.split('_')[1]
# get itemized information for the DataFrame
exp = infile.replace('IC_','')
exp = exp.replace('.p','')
basin, year, season = exp.split('_')

seg2_list = flux_fun.ic_seg2_dict[source]
    
aa = pd.read_pickle(indir + infile)

this_aa = aa.loc[:,seg2_list]
this_V = V[seg2_list]
net_V = this_V.sum()
this_rv = this_V / net_V # relative volume

this_net_aa = this_aa.copy()

for sn in this_V.index:
    VV = this_V[sn]
    this_net_aa.loc[:,sn] = this_net_aa.loc[:,sn] * VV

# make a Series of mean concentration in the volume
mean_c = this_net_aa.sum(axis=1) / net_V

# find e-folding time
td = mean_c.index.to_numpy()
mc = mean_c.to_numpy()
ind_ef = np.argwhere(mc < 1/np.e)[0]
tres = td[ind_ef] # residence time in days

# make a mask of segment names
acol = list(aa.columns)
not_seg2_list = [item for item in acol if item not in seg2_list]
print(' %10s: Tres = %0.2f days' % (exp, tres))

# residence time figure
plt.close('all')
fs = 16
plt.rc('font', size=fs)

fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)

for s in seg2_list:
    if False:#s == 'G4_s':
        pass
        #aa[s].plot(ax=ax, legend=False, ls='-', c='r', lw=50*this_rv[s], alpha=.5)
    else:
        #aa[s].plot(ax=ax, legend=False, ls='-', c='gold', lw=50*this_rv[s], alpha=.5)
        aa[s].plot(ax=ax, ls='-', lw=50*this_rv[s], alpha=.5, legend=True)
# for s in not_seg2_list:
#     if s == 'J4_f':
#         aa[s].plot(ax=ax, legend=False, ls='-', c='g', alpha=.5)
#     elif s == 'J4_s':
#         aa[s].plot(ax=ax, legend=False, ls='-', c='b', alpha=.5)
#     else:
#         aa[s].plot(ax=ax, legend=False, ls='-', c='c', alpha=.5)


ax.plot(tres,1/np.e,'*k', ms=30)

fig.tight_layout()
#fig.savefig(outdir + 'junk.png')
plt.show()
plt.rcdefaults()


