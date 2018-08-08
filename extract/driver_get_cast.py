#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:22:34 2018

@author: pm7
"""

import pandas as pd

import cast_fun as cfun
from importlib import reload
reload(cfun)

# set defaults
gridname = 'cas4'
tag = 'v2'
ex_name = 'lo6biom'

testing = False

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

# +++ load ecology CTD cast data +++

dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# add Canadian data

dir1 = Ldir['parent'] + 'ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')

sta_df = pd.concat((sta_df, sta_df_ca))

year = 2017
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')

Casts = pd.concat((Casts, Casts_ca))


# trim the station list as desired
sta_list = [sta for sta in sta_df.index]

for station in sta_list:
    
    casts = Casts[Casts['Station'] == station]
    casts = casts.set_index('Date')
    casts = casts.loc[:,['Salinity', 'Temperature','Z']] # keep only selected columns
    # identify a single cast by its date
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    
    lon_str = str(sta_df.loc[station,'Longitude'])
    lat_str = str(sta_df.loc[station,'Latitude'])

    ds_list = []
    for cdate in castdates:
        ds_list.append(cdate.strftime('%Y.%m.%d'))
        if testing:
            print(station + ' ' + cdate.strftime('%Y.%m.%d') + ' ' + lon_str + ' ' + lat_str)
        
    if not testing:
        for date_string in ds_list:   
            cfun.get_cast(gridname, tag, ex_name, date_string, station, lon_str, lat_str)