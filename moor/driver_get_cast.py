#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:22:34 2018

@author: pm7
"""

import cast_fun as cfun
from importlib import reload
reload(cfun)

# set defaults
gridname = 'cas3'
tag = 'v0'
ex_name = 'lo6m'
sta_name = 'HCB010'
lon_str = '-122.8200'
lat_str = '47.6670'

ds_list = [
        '2017.01.12',
        '2017.02.03',
        '2017.03.09',
        '2017.04.06',
        '2017.05.18',
        '2017.06.08',
        '2017.07.06',
        '2017.08.11',
        '2017.09.06',
        '2017.10.05',
        '2017.11.02',
        '2017.12.13',
        ]

for date_string in ds_list:   
    cfun.get_cast(gridname, tag, ex_name, date_string, sta_name, lon_str, lat_str)