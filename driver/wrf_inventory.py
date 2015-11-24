"""
This code automatically populates a pandas dataframe with information about
the WRF files in our archive.

Designed to be run only on fjord.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd

if Ldir['parent'] == '/Users/PM5/Documents/':
    f_dir0 = Ldir['parent'] + 'LiveOcean_files/wrf/'
elif Ldir['parent'] == '/data1/parker/': 
    f_dir0 = '/pmraid3/darr/tstwrf/tmpwrf/'

# directory names are like 2012100700

# expected directories
from datetime import datetime, timedelta
fe_list = []
dt0 = datetime(2012,10,7)
dt = dt0
while datetime.now() - dt > timedelta(1):
    fe_list.append(dt.strftime('%Y%m%d') + '00')
    dt = dt + timedelta(1)
    
clist = ['dir_exists','N_d2','N_d3']
f_df = pd.DataFrame(index = fe_list, columns = clist)

# what forecasts exist in the forcing directory
f_list = []
for item in os.listdir(f_dir0):
    if item[-2:] == '00' and len(item) == 10:
        f_list.append(item)
        
f_df['dir_exists'] = 'no'

for item in f_list:
    if item in fe_list:
        f_df.ix[item,'dir_exists'] = 'YES'
        
df_yes = f_df[f_df['dir_exists']=='YES']
df_no = f_df[f_df['dir_exists']=='no']

for item in df_yes.index:
    fl = os.listdir(f_dir0 + item)
    # how many CONTINUOUS forecasts have been saved successfully
    # file names are like wrfout.ocean_d2.2012100700.f00.0000
    flag = 1
    nd = 0
    while flag == 1:
        sd = ('00' + str(nd))[-2:]
        if ('wrfout.ocean_d2.' + item + '.f' + sd + '.0000') in fl:
            nd += 1
        else:
            flag = 0
    df_yes.ix[item]['N_d2'] = nd
       
    flag = 1
    nd = 0
    while flag == 1:
        sd = ('00' + str(nd))[-2:]
        if ('wrfout.ocean_d3.' + item + '.f' + sd + '.0000') in fl:
            nd += 1
        else:
            flag = 0
    df_yes.ix[item]['N_d3'] = nd

df_yes_but =  df_yes[(df_yes['N_d2'] < 25) | (df_yes['N_d3'] < 25)]

# make sure that the dates are in order
df_yes_but = df_yes_but.sort_index()
df_no = df_no.sort_index()

# print to the screen
pd.set_option('display.max_rows', 1000)
print('MISSING DIRECTORIES')
print(df_no)
print('\nMISSING HOURS')
print(df_yes_but)






