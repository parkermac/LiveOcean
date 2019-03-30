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

# make vector and array times in days from start of the year
dt = []
for tt in ot:
    dt.append(Lfun.modtime_to_datetime(tt))
    
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
