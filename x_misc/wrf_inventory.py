"""
This code automatically populates a pandas dataframe with information about
the WRF files in our archive.

"""

# setup
import os; import sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd

if Ldir['lo_env'] == 'pm_mac':
    f_dir0 = Ldir['parent'] + 'LiveOcean_data/wrf/'
else: 
    f_dir0 = '/pmr2/darr/wrf_crons/wrfout/'

# directory names are like 2012100700

# expected directories
from datetime import datetime, timedelta
fe_list = []
dt0 = datetime(2012,10,7)
dt = dt0
while datetime.now() - dt > timedelta(1):
    fe_list.append(dt.strftime('%Y%m%d') + '00')
    dt = dt + timedelta(1)
    
clist = ['dir_exists','N_d2','N_d3','N_d4']
f_df = pd.DataFrame(index = fe_list, columns = clist)

# what forecasts exist in the forcing directory
f_list = []
LL = os.listdir(f_dir0)
LL.sort()
for item in LL:
    if item[-2:] == '00' and len(item) == 10:
        f_list.append(item)
        
f_df['dir_exists'] = 'no'

for item in f_list:
    if item in fe_list:
        f_df.loc[item,'dir_exists'] = 'YES'
        
df_yes = f_df[f_df['dir_exists']=='YES']
df_no = f_df[f_df['dir_exists']=='no']

for item in df_yes.index:
    fl = os.listdir(f_dir0 + item)
    fl.sort()
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
    df_yes.loc[item]['N_d2'] = nd
       
    flag = 1
    nd = 0
    while flag == 1:
        sd = ('00' + str(nd))[-2:]
        if ('wrfout.ocean_d3.' + item + '.f' + sd + '.0000') in fl:
            nd += 1
        else:
            flag = 0
    df_yes.loc[item]['N_d3'] = nd
    
    flag = 1
    nd = 0
    while flag == 1:
        sd = ('00' + str(nd))[-2:]
        if ('wrfout.ocean_d4.' + item + '.f' + sd + '.0000') in fl:
            nd += 1
        else:
            flag = 0
    df_yes.loc[item]['N_d4'] = nd
    

df_yes_but =  df_yes[(df_yes['N_d2'] < 25) | (df_yes['N_d3'] < 25) | (df_yes['N_d4'] < 25)]

# make sure that the dates are in order
df_yes_but = df_yes_but.sort_index()
df_no = df_no.sort_index()

# print to the screen
pd.set_option('display.max_rows', 1000)
print('MISSING DIRECTORIES')
print(df_no)
print('\nMISSING HOURS')
print(df_yes_but)






