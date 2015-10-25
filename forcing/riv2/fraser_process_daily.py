# -*- coding: utf-8 -*-
"""
This concatenates all the daily Fraser files pulled in from gmail.
It makes separate files for flow and temperature.

Performance: 11 sec to process 600 days of files.

Here is example header info:
    
Environment Canada, Real-time hydrometric data, Copyright 2014
Report generated at: 2014/02/12 11:19:00
Report type: raw data values

Station ID, Station Name
08MF005,FRASER RIVER AT HOPE

Variable ID, Variable Name
47,Discharge provisional  (m3/s)
5,Water temperature (°C)

"""

import os; import sys; pth = os.path.abspath('../../alpha')
if pth not in sys.path: sys.path.append(pth)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart()

import pandas as pd

indir = Ldir['data'] + 'rivers/data_fraser_raw/'
outdir = Ldir['data'] + 'rivers/data_processed/'

# create a list of files
flist = os.listdir(indir)
fl = []
for fn in flist:
    if 'Parker' in fn:
        fl.append(fn)

qdf_list = []
tdf_list = []
         
for fn in fl:
    try:
        rdf = pd.read_csv(indir + fn, skiprows=11)
        # Note: hard coding the skiprows is brittle, but you never know
        # what might change about the formatting, and this is one way
        # to potentially be alerted to changes.
    except: # raises CParserError but this is not recognized
        print('Warning: could not read ' + fn)
        continue    
    cols = rdf.columns
    new_cols = []
    for cc in cols:
        new_cols.append(cc.strip().lower())
    rdf.columns = new_cols
    rdf = rdf.set_index('date')
    
    # pull the flow into its own Series
    qdf = rdf[rdf['parameter identification']==47]
    qdf = qdf['value']
    qdf.index = qdf.index.to_datetime()
    
    # pull the temperature into its own Series
    tdf = rdf[rdf['parameter identification']==5]
    tdf = tdf['value']
    tdf.index = tdf.index.to_datetime()
    
    qdf_list.append(qdf)
    tdf_list.append(tdf)

# pull all series together   
qdf = pd.concat(qdf_list)
tdf = pd.concat(tdf_list)

# make daily averages (data stamp will be noon)
qdf = qdf.resample('D', how='mean', label='right', loffset='-12h')
tdf = tdf.resample('D', how='mean', label='right', loffset='-12h')

# remove bad values
qdf = qdf[qdf.notnull()]
tdf = tdf[tdf.notnull()]

# save to files
qdf.to_pickle(outdir + 'fraser_flow.p')
tdf.to_pickle(outdir + 'fraser_temperature.p')

