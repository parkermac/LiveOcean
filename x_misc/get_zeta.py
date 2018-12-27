#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract zeta time series from regions of the hycom forcing.  The goal is to develop background information in support of an improved open boundary condition in the Northern Strait of Georgia.

Note that Days 200-250 are a fairly stable time to explore a solution (2017.07.20 to 2017.09.08).  This is also a time period where DT was generally long (0 or 1 blowups).

"""

from datetime import datetime, timedelta
import argparse
import pandas as pd
import pickle
import numpy as np

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun

# Command line arguments
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2018.11.26')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2018.12.25')
args = parser.parse_args()

# use some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)
dt0 = datetime.strptime(args.date_string0, '%Y.%m.%d')
dt1 = datetime.strptime(args.date_string1, '%Y.%m.%d')

# make a list of model date strings
mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, '%Y.%m.%d'))
    mdt = mdt + timedelta(days=1)

# extractions
zdf = pd.DataFrame(columns=['za', 'zb'])

print('== looping over dates ==')
# loop over dates
for mds in mds_list:

    indir0 = Ldir['LOo'] + Ldir['gtag'] + '/f' + mds + '/'

    # specify which ocn forcing to look in
    get_data = True
    if os.path.isdir(indir0 + 'ocn3'):
        ocn = 'ocn3'
        print(' - ' + mds + ': using ocn3')
    elif os.path.isdir(indir0 + 'ocn2'):
        ocn = 'ocn2'
        print(' - ' + mds + ': using ocn2')
    else:
        print(' - ' + mds + ': no ocn directory')
        get_data = False

    if get_data:
        
        # get coordinates
        indir = indir0 + ocn + '/Data/'
        coord_dict = pickle.load(open(indir + 'coord_dict.p', 'rb'))
        lon = coord_dict['lon']
        lat = coord_dict['lat']

        # get fields
        # note that when we use the "xfh" fields there should be no nans
        xfh = pickle.load(open(indir + '/xfh' + mds + '.p', 'rb'))
        zeta = xfh['ssh']

        # get indices around certain regions, and their data
        #
        # (a) Mouth of Strait of Juan de Fuca
        i0, i1, fr = zfun.get_interpolant(np.array([-124.7, -124.5]), lon, extrap_nan=False)
        j0, j1, fr = zfun.get_interpolant(np.array([48.4, 48.6]), lat, extrap_nan=False)
        zeta_a = zeta[j0[0]:j1[1]+1, i0[0]:i1[1]+1]
        #
        # (b) Northern Strait of Georgia
        i0, i1, fr = zfun.get_interpolant(np.array([-125.5, -124.5]), lon, extrap_nan=False)
        j0, j1, fr = zfun.get_interpolant(np.array([50.0, 50.3]), lat, extrap_nan=False)
        zeta_b = zeta[j0[0]:j1[1]+1, i0[0]:i1[1]+1]
        
        mdt = datetime.strptime(mds,'%Y.%m.%d')
        zdf.loc[mdt, 'za'] = zeta_a.mean()
        zdf.loc[mdt, 'zb'] = zeta_b.mean()
        
# save output
outdir = Ldir['LOo'] + 'misc/'
Lfun.make_dir(outdir)
zdf.to_pickle(outdir + 'zeta_df.p')
        

