"""
Plot exchange flow vs. ds/dx at a sepcific section as a scatterplot
for all years and seasons.

"""

# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
gridname = 'cas6'; tag = 'v3'
Ldir = Lfun.Lstart(gridname, tag)


import tef_fun
import flux_fun

# select input directory
indir0 = Ldir['LOo'] + 'tef/'

qe_list = []
ds_list = []
for year in [2017, 2018, 2019]:
    year_str = str(year)
    indir = indir0 + 'cas6_v3_lo8b_'+year_str+'.01.01_'+year_str+'.12.31/flux/'
    for season in ['winter', 'spring', 'summer', 'fall']:
        fn = indir + 'two_layer_' + season + '.p'
        df = pickle.load(open(fn,'rb'))

        q_s = df.loc['ai2','q_s']/1e3
        q_f = df.loc['ai2','q_f']/1e3
        s_s1 = df.loc['ai1','s_s']
        s_f1 = df.loc['ai1','s_f']
        s_s3 = df.loc['ai3','s_s']
        s_f3 = df.loc['ai3','s_f']
        
        qe_list.append( (np.abs(q_s)+np.abs(q_s))/2 )
        ds_list.append( (s_s1+s_f1)/2 - (s_s3+s_f3)/2 )
            
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(ds_list, qe_list, 'ob')

plt.show()
    

    