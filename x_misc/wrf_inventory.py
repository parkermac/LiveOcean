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
from datetime import datetime, timedelta

if Ldir['lo_env'] == 'pm_mac':
    f_dir0 = Ldir['data'] + 'wrf/'
else: 
    f_dir0 = '/pmr2/darr/wrf_crons/wrfout/'


#expected dates
dte = pd.date_range(start=datetime(2012,10,7), end=datetime.now())
    
clist = ['dir_name', 'dir_exists', 'd2', 'd3', 'd4']

f_df = pd.DataFrame(0, index = dte, columns = clist)

for dt in dte:
    # directory names are like 2012100700
    f_df.loc[dt, 'dir_name'] = dt.strftime('%Y%m%d' + '00')

# what 00 forecasts exist in the forcing directory
f_list = os.listdir(f_dir0)
f_list = [item for item in f_list if (item[-2:] == '00' and len(item) == 10)]
f_list.sort()

# list of datetimes of forecasts that exist
fdt_list = [datetime.strptime(item[:8],'%Y%m%d') for item in f_list]
        
# note which dates we have forecast directories for them, and which do not
f_df.loc[:,'dir_exists'] = 'no'
f_df.loc[fdt_list, 'dir_exists'] = 'yes'

for dt in fdt_list:
    # count up the number of each type of file on each day
    
    dname = f_df.loc[dt, 'dir_name']
    
    fl = os.listdir(f_dir0 + dname)
    fl.sort()
    
    n2 = len([item for item in fl if 'd2' in item])
    n3 = len([item for item in fl if 'd3' in item])
    n4 = len([item for item in fl if 'd4' in item])
    
    f_df.loc[dt, 'd2'] = n2
    f_df.loc[dt, 'd3'] = n3
    f_df.loc[dt, 'd4'] = n4
    
f_df_no = f_df[f_df['dir_exists']=='no']
f_df_yes = f_df[f_df['dir_exists']=='yes']

print(50*'=')
print('Number of days with missing forecasts = %d' % (len(f_df_no)))
print('\n' + 50*'=')

print('*** Monthly Statistics ***')
# print out monthly statistics
pd.set_option('display.max_rows', 1000)
fm_df_yes = f_df_yes.resample('M').mean()
print(fm_df_yes)

# NOTE: you could print out a range of daily information with a command like
# f_df_yes.loc[pd.date_range(datetime(2018,1,1),datetime(2018,1,10)),:]




