"""
Process historical Fraser Data.

(1) Daily__Oct-18-2015_10_41_57PM__ddf.csv

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

++++++++++++++++++++++++++++++++++++++++++++++

(2) 08mf005_prelim_2013_2014.csv

The second file of 2013-2014 data that I got from Environment Canada
looks like:
    
Unapproved data - subject to revision,,,,,,,,,,
Discharge.Working@08MF005,,08MF005,Filter=None,,m^3/s,,,,,
FRASER RIVER AT HOPE,,,,,,,,,,
Day,Mean,Grade,,,,,,,,
1/1/13,876,,,,,,,,,
1/2/13,876,,,,,,,,,
1/3/13,856,,,,,,,,,
1/4/13,844,,,,,,,,,
1/5/13,863,,,,,,,,,

"""

import os; import sys; pth = os.path.abspath('../../alpha')
if pth not in sys.path: sys.path.append(pth)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart()

import pandas as pd

indir = Ldir['data'] + 'rivers/data_fraser_historical/'
outdir = Ldir['data'] + 'rivers/data_processed/'

# File (1)

fn = 'Daily__Oct-18-2015_10_41_57PM__ddf.csv'
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

# File (2)

fn = '08mf005_prelim_2013_2014.csv'
rdf = pd.read_csv(indir + fn, skiprows=3)
cols = rdf.columns
new_cols = []
for cc in cols:
    new_cols.append(cc.strip().lower())
rdf.columns = new_cols
rdf = rdf.set_index('day')
qdf = rdf['mean']
qdf.index = qdf.index.to_datetime()
# make daily averages (data stamp will be noon)
qdf = qdf.resample('D', how='mean', label='right', loffset='-12h')
# remove bad values
qdf = qdf[qdf.notnull()]
# save to files
qdf.to_pickle(outdir + 'fraser_flow_historical_2.p')

# File (3) Temperature

# This data is still in need of cleaning, so I will ignore it for now.

if False:
    fn = '08mf005_prelim_2008_2015_wtemp.csv'
    rdf = pd.read_csv(indir + fn, skiprows=3)
    cols = rdf.columns
    new_cols = []
    for cc in cols:
        new_cols.append(cc.strip().lower())
    rdf.columns = new_cols
    rdf = rdf.set_index('day')
    qdf = rdf['mean']
    qdf.index = qdf.index.to_datetime()
    # make daily averages (data stamp will be noon)
    #qdf = qdf.resample('D', how='mean', label='right', loffset='-12h')
    # remove bad values
    qdf = qdf[qdf.notnull()]
    # save to files
    qdf.to_pickle(outdir + 'fraser_temperature_historical.p')
