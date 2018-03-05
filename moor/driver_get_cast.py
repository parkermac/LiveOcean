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
gridname = 'cas3'
tag = 'v0'
ex_name = 'lo6m'

testing = False

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

# where the data is
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
fn = dir0 + 'ParkerMacCready2017CTDDataFeb2018.xlsx'

# load station location and depth info
sta_info_fn = dir0 + 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx'
sta_df = pd.read_excel(sta_info_fn)
sta_df = sta_df.set_index('Station')
# get locations in decimal degrees
for sta in sta_df.index:
    lat_str = sta_df.loc[sta, 'Lat_NAD83 (deg / dec_min)']
    lat_deg = float(lat_str.split()[0]) + float(lat_str.split()[1])/60
    sta_df.loc[sta,'Latitude'] = lat_deg
    #
    lon_str = sta_df.loc[sta, 'Long_NAD83 (deg / dec_min)']
    lon_deg = float(lon_str.split()[0]) + float(lon_str.split()[1])/60
    sta_df.loc[sta,'Longitude'] = -lon_deg    
sta_df.pop('Lat_NAD83 (deg / dec_min)')
sta_df.pop('Long_NAD83 (deg / dec_min)')

# read in the data (all stations, all casts)
all_casts = pd.read_excel(fn, sheet_name='2017Provisional_CTDResults',
                          parse_dates = ['Date'])
# trim the station list uf desired
if testing:
    sta_list = [sta for sta in sta_df.index if 'HCB' in sta]
else:
    sta_list = [sta for sta in sta_df.index]

for station in sta_list:
    print('Extracting: ' + station)           
    casts = all_casts[all_casts['Station'] == station]   
    casts = casts.set_index('Date')    
    
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    
    title = station + ': ' + sta_df.loc[station,'Descrip']
    Max_z = -float(sta_df.loc[station, 'Max_Depth'])
    
    lon_str = str(sta_df.loc[station,'Longitude'])
    lat_str = str(sta_df.loc[station,'Latitude'])

    ds_list = []
    for cdate in castdates:
        ds_list.append(cdate.strftime('%Y.%m.%d'))
        
    if testing:
        ds_list = ['2017.03.30']

    for date_string in ds_list:   
        cfun.get_cast(gridname, tag, ex_name, date_string, station, lon_str, lat_str)