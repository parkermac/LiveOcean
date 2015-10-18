"""
Process historical Fraser Data.

The start of the file looks like:
    
Daily Discharge and Daily Water Level
  ID,PARAM,Date,Value,SYM
08MF005,1,2000/01/01,1080,
08MF005,1,2000/01/02,1080,
08MF005,1,2000/01/03,1070,
08MF005,1,2000/01/04,1060,
08MF005,1,2000/01/05,1060,
08MF005,1,2000/01/06,1010,
08MF005,1,2000/01/07,981,
08MF005,1,2000/01/08,987,
08MF005,1,2000/01/09,1000,
08MF005,1,2000/01/10,1030,
08MF005,1,2000/01/11,1010,
08MF005,1,2000/01/12,973,
08MF005,1,2000/01/13,939,A
08MF005,1,2000/01/14,909,

I'm not sure what "A" means.
"""

import os; import sys; pth = os.path.abspath('../../alpha')
if pth not in sys.path: sys.path.append(pth)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart()

import pandas as pd

indir = Ldir['data'] + 'rivers/data_fraser_historical/'
outdir = Ldir['data'] + 'rivers/data_processed/'

fn = 'Daily__Oct-18-2015_06_30_13PM__ddf.csv'

rdf = pd.read_csv(indir + fn, skiprows=1)
cols = rdf.columns
new_cols = []
for cc in cols:
    new_cols.append(cc.strip().lower())
rdf.columns = new_cols
rdf = rdf.set_index('date')

# pull the flow into its own Series
qdf = rdf[rdf['param']==1]
qdf = qdf['value']
qdf.index = qdf.index.to_datetime()

# make daily averages (data stamp will be noon)
qdf = qdf.resample('D', how='mean', label='right', loffset='-12h')

# remove bad values
qdf = qdf[qdf.notnull()]

# save to files
qdf.to_pickle(outdir + 'fraser_flow_historical.p')
