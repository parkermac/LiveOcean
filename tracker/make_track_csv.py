# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:43:43 2016

@author: PM5
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pandas as pd

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
inname = m_dict[my_npt]
  
import pickle  # python 3
P, G, S, PLdir = pickle.load( open( indir + inname, 'rb' ) )

# make sure the output directory exists
outdir = Ldir['LOo'] + 'tracks/'
outname = inname[:-2]
odir = outdir + outname + '/'
Lfun.make_dir(odir)

for vn in P.keys():
    a = pd.DataFrame(P[vn], index=P['ot'])
    a.index.name = 'seconds_since_1_1_1970'
    a.to_csv(odir + vn)