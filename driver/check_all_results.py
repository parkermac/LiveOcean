"""
This code automatically populates a pandas dataframe with information about
the status of LiveOcean.  For each forecast day, what forcing files have been 
created?  Has the day run successfully? Etc.
"""

# get command line arguments if any
import argparse
parser = argparse.ArgumentParser()
# optional arguments
parser.add_argument("-g", "--gridname", type=str, default='cascadia1', help="cascadia1, etc.")
parser.add_argument("-t", "--tag", type=str, default='base', help="base, etc.")
parser.add_argument("-x", "--ex_name", type=str, default='lo1', help="e.g. lo1")
args = parser.parse_args()
# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

import pandas as pd

# what forecasts exist in the forcing directory
f_dir0 = Ldir['LOo'] + Ldir['gtag'] + '/'
f_dir0_list = []
for item in os.listdir(f_dir0):
    if item[0] == 'f' and len(item) == 11:
        f_dir0_list.append(item)
        
force_dict = {'atm':['lwrad_down.nc','Pair.nc','Qair.nc','rain.nc','swrad.nc',
            'Tair.nc','Uwind.nc','Vwind.nc'],
        'ocn':['ocean_bry.nc','ocean_clm.nc','ocean_ini.nc'],
        'riv':['rivers.nc',],
        'tide':['tides.nc',]}
    
clist = ['atm','ocn','riv','tide','rundef','his']
f_df = pd.DataFrame(index = f_dir0_list, columns = clist)

for which_forecast in f_df.index:
    for which_force in force_dict.keys():
        force_dir = f_dir0 + which_forecast + '/' + which_force
        try:
            lll = os.listdir(force_dir)
            nc_list = force_dict[which_force]
            if set(nc_list).issubset(set(lll)):
                f_df.ix[which_forecast, which_force] = 'YES'
            else:
                f_df.ix[which_forecast, which_force] = 'no'
        except:
            # assume the directory is missing
            f_df.ix[which_forecast, which_force] = '--'

# what .in files exist
r_dir0 = Ldir['roms'] + 'output/' + Ldir['gtag'] + '/'
r_dir0_list = []
for item in os.listdir(r_dir0):
    if item[0] == 'f' and len(item) == 11:
        f_string = item
        if 'liveocean.in' in os.listdir(r_dir0 + f_string):
            f_df.ix[f_string,'rundef'] = 'YES'
        
# what forecasts have been run successfully
out_list = []
for ii in range(2,26):
    ncpad = '0000' + str(ii)
    ncpad = ncpad[-4:]
    hisname = 'ocean_his_' + ncpad + '.nc'
    out_list.append(hisname)
ro_dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
for item in os.listdir(ro_dir0):
    if item[0] == 'f' and len(item) == 11:
        if set(out_list).issubset(set(os.listdir(ro_dir0 + item))):
                    f_df.ix[item, 'his'] = 'YES'
                    
# eventually we would also like to check what has been pushed to azure

# mark missing things
f_df[f_df.isnull()] = '--'

# make sure that the dates are in order
f_df = f_df.sort_index()

# print to the screen
pd.set_option('display.max_rows', 1000)
print f_df

if False:
    # and save in an html file
    from datetime import datetime
    fn = open(Ldir['LOo'] + 'forecast_state_' +
        datetime.now().strftime('%Y.%m.%d') + '.html','w')
    fn.write(f_df.to_html())
    fn.close()





