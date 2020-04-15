#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to make cast-like extractions from a LiveOcean run, at times
and places that match an observational dataset,
specified by -ds [ecology, woac, ...] and a year.
"""
import os; import sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import pandas as pd
import cast_fun as cfun
from importlib import reload
reload(cfun)

testing = False

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas6')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v3')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo8b')
parser.add_argument('-y', '--year', nargs='?', type=int, default=2017)
parser.add_argument('-ds', '--data_source', nargs='?', type=str, default='ecology')
args = parser.parse_args()

gridname = args.gridname
tag = args.tag
ex_name = args.ex_name
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name

year = args.year
data_source = args.data_source

if data_source == 'ecology_canada':
    # +++ load ecology CTD cast data +++
    dir0 = Ldir['parent'] + 'ptools_data/ecology/'
    # load processed station info and data
    sta_df = pd.read_pickle(dir0 + 'sta_df.p')
    # add Canadian data
    dir1 = Ldir['parent'] + 'ptools_data/canada/'
    # load processed station info and data
    sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
    # combine
    sta_df = pd.concat((sta_df, sta_df_ca), sort=False)
    # and get cast data
    Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
    Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
    #
    Casts = pd.concat((Casts, Casts_ca), sort=False)
if data_source == 'ecology':
    # +++ load ecology CTD cast data +++
    dir0 = Ldir['parent'] + 'ptools_data/ecology/'
    # load processed station info and data
    sta_df = pd.read_pickle(dir0 + 'sta_df.p')
    # and get cast data
    Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
elif data_source == 'woac':
    dir0 = Ldir['parent'] + 'ptools_data/' + data_source + '/'
    # load processed station info and data
    sta_df = pd.read_pickle(dir0 + 'sta_df.p')
    # and get cast data
    Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
    
# trim the station list as desired
sta_list = [sta for sta in sta_df.index if 'SOG' in sta]

if 'ecology' in data_source:
    for station in sta_list:
        casts = Casts[Casts['Station'] == station]
        casts = casts.set_index('Date')
        casts = casts.loc[:,['Salinity', 'Temperature','Z']] # keep only selected columns
        # identify a single cast by its date
        alldates = casts.index
        castdates = alldates.unique() # a short list of unique dates (1 per cast)
        #
        lon_str = str(sta_df.loc[station,'Longitude'])
        lat_str = str(sta_df.loc[station,'Latitude'])
        #
        ds_list = []
        for cdate in castdates:
            ds_list.append(cdate.strftime('%Y.%m.%d'))
            if testing:
                print(station + ' ' + cdate.strftime('%Y.%m.%d') + ' ' + lon_str + ' ' + lat_str)
        #
        if not testing:
            for date_string in ds_list:   
                cfun.get_cast(gridname, tag, ex_name, date_string, station, lon_str, lat_str)
elif data_source == 'woac':
    for cast in sta_df.index:
        station = str(sta_df.loc[cast,'Station'])
        date_str = sta_df.loc[cast,'Datetime'].strftime('%Y.%m.%d')
        lon_str = str(sta_df.loc[cast,'Longitude'])
        lat_str = str(sta_df.loc[cast,'Latitude'])
        cfun.get_cast(gridname, tag, ex_name, date_str, 'WOAC'+station, lon_str, lat_str)
        print('%s %s %s %s' % (date_str, 'WOAC'+station, lon_str, lat_str))
