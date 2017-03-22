#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:25:59 2017

@author: PM5

Code to test getting small collections of hycom1 files.
"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname='cas1', tag='base')

from datetime import datetime, timedelta

# note: there is a big gap between:
# h2012.09.03.p and
# h2012.12.02.p

Ldir['run_type'] = 'backfill'
#Ldir['date_string'] = '2012.09.19'
Ldir['date_string'] = '2012.12.01'

in_dir = Ldir['data'] + 'hycom1/'

# get a list of all available times
h_list0 = os.listdir(in_dir)
h_list = [item for item in h_list0 if item[0] == 'h']
# then find the index of the start of the current day
# but if it is missing search for the most recent past one that exists
keep_looking = True
dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
it0 = None
counter = 0
maxcount = 100 # this handles the biggest gap we have
while keep_looking and counter < maxcount:
    dt_next = dt_now - timedelta(days=counter)
    dts_next = datetime.strftime(dt_next, '%Y.%m.%d')
    try:
        it0 = h_list.index('h' + dts_next + '.p')
        keep_looking = False
        if counter > 0:
            print('Warning: Needed %d iterations' % (counter))
    except ValueError:
        counter += 1       
# print the list of files
if it0 == None:
    print('ERROR: no valid files found at nearby times')    
else:        
    it_list = range(it0-2, it0+4)
    for it in it_list:
        fn = h_list[it]
        print(fn) # debugging
        
        
