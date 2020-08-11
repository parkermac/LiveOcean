"""
Plot profiles of salt or density vs. depth at each section.
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
indir = indir0 + 'profiles/'

sect_list = ['ai1','ai3']

plt.close('all')
fig = plt.figure(figsize=(10,10))
ax_list = []
for mo in range(1,13,1):
    ax_list.append(fig.add_subplot(3,4,mo))


df_dict = {}
for sect in sect_list:
    df_dict[sect] = pd.read_pickle(indir + sect + '.p')

z = df_dict[sect_list[0]].index.to_numpy()
for mo in range(1,13,1):
    for sect in sect_list:
        this_s = df_dict[sect].loc[:,mo].to_numpy()
        ax_list[mo-1].plot(this_s, z)
        
for mo in range(1,13,1):
    ax_list[mo-1].set_xlim(28,33)
        
plt.show()