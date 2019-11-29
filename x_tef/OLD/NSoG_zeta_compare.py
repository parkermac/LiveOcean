"""
Plot SSH (zeta) from a TEF extraction, compared with the value
found in the bry forcing.  Hard-wired here to look at TEF section
sog5 and the northern boundary.

"""
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

import numpy as np
import matplotlib.pyplot as plt
import pickle
from datetime import datetime, timedelta
import pandas as pd

import tef_fun
# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# choose input and organize output
Ldir = Lfun.Lstart()
indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to process
item = Lfun.choose_item(indir0)
indir0 = indir0 + item + '/'
indir = indir0 + 'bulk/'

in_fn = indir + 'sog5.p'

bulk = pickle.load(open(in_fn, 'rb'))
ot = bulk['ot']
ssh = bulk['ssh_lp'] # average SSH across the section, low-passed
qnet = bulk['qnet_lp'] # net transport through the section, low-passed

# make vector and array times in days from start of the year
dt = []
for tt in ot:
    dt.append(Lfun.modtime_to_datetime(tt) - timedelta(-0.5))
    
# find the correspinding bry extraction
a = indir0.split('/')[-2]
aa = a.split('_')
gridname = aa[0]
tag = aa[1]
ex_name = aa[2]
ds0 = aa[3]
ds1 = aa[4]

bry_fn = Ldir['LOo'] + 'misc/bry_' + gridname + '_' + tag + '_' + ds0 + '_' + ds1 + '/zeta_north.p'

bry = pd.read_pickle(bry_fn)

bry_df = pd.DataFrame(index = bry.index, columns = ['zeta_sog5', 'zeta_bry', 'dzeta', 'qnet'])
bry_df['zeta_bry'] = bry
bry_df.loc[dt,'zeta_sog5'] = ssh
bry_df.loc[dt,'qnet'] = qnet
bry_df['dzeta'] = bry_df['zeta_sog5'] - bry_df['zeta_bry']
tstr = gridname + '_' + tag

fig = plt.figure(figsize = (14,14))

ax = fig.add_subplot(221)
bry_df.plot(y=['zeta_sog5','zeta_bry'], ax=ax, title=tstr,
    xlim=(datetime(2017,1,1), datetime(2017,12,31)), ylim=(0,0.7), grid=True)
    
ax = fig.add_subplot(222)
bry_df.plot(x='dzeta', y='qnet', ax=ax, grid=True, style='*r')
ax.text(.1, .9, 'zeta sog5 - bry', transform=ax.transAxes)
ax.axhline()
ax.axvline()
    
ax = fig.add_subplot(223)
bry_df.plot(x='zeta_sog5', y='qnet', ax=ax, grid=True, style='*b')
ax.text(.1, .9, 'zeta sog5', transform=ax.transAxes)
ax.axhline()
ax.axvline()
    
ax = fig.add_subplot(224)
bry_df.plot(x='zeta_bry', y='qnet', ax=ax, grid=True, style='*g')
ax.text(.1, .9, 'zeta bry', transform=ax.transAxes)
ax.axhline()
ax.axvline()

plt.show()
