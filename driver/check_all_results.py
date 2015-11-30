"""
This code automatically populates a pandas dataframe with information about
the status of LiveOcean.  For each forecast day, what forcing files have been
created?  Has the day run successfully? Etc.
"""

# setup
import os
import sys
import argparse
import pandas as pd
from datetime import datetime, timedelta

alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

# get command line arguments, if any
parser = argparse.ArgumentParser()
# optional arguments
parser.add_argument("-g", "--gridname", type=str, default='cascadia1')
parser.add_argument("-t", "--tag", type=str, default='base')
parser.add_argument("-x", "--ex_name", type=str, default='lo1')
parser.add_argument("-nd", "--num_days", type=int, default=30)
args = parser.parse_args()

Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

# which forecast days exist in the forcing directory
f_dir0 = Ldir['LOo'] + Ldir['gtag'] + '/'
f_dir0_list = []
for item in os.listdir(f_dir0):
    if item[0] == 'f' and len(item) == 11:
        f_dir0_list.append(item)

force_dict = {'atm': ['lwrad_down.nc', 'Pair.nc', 'Qair.nc', 'rain.nc',
                      'swrad.nc', 'Tair.nc', 'Uwind.nc', 'Vwind.nc'],
              'ocn': ['ocean_bry.nc', 'ocean_clm.nc', 'ocean_ini.nc'],
              'riv': ['rivers.nc', ],
              'tide': ['tides.nc', ]}

clist = ['atm', 'ocn', 'riv', 'tide', 'dot_in', 'his', 'lp', 'azu']

f_df = pd.DataFrame(index=f_dir0_list, columns=clist)

for which_forecast in f_df.index:
    for which_force in force_dict.keys():
        force_dir = f_dir0 + which_forecast + '/' + which_force + '/'
        try:
            lll = os.listdir(force_dir)
            nc_list = force_dict[which_force]
            if set(nc_list).issubset(set(lll)):
                f_df.ix[which_forecast, which_force] = 'YES'
                if which_force in ['atm', 'ocn']:
                    try:
                        ps = Lfun.csv_to_dict(force_dir + 'Info/process_status.csv')
                        dt0 = datetime.strptime(ps['start_time'], '%d-%b-%Y %H:%M:%S')
                        dt1 = datetime.strptime(ps['end_time'], '%d-%b-%Y %H:%M:%S')
                        vdt0 = datetime.strptime(ps['var_start_time'], '%d-%b-%Y')
                        vdt1 = datetime.strptime(ps['var_end_time'], '%d-%b-%Y')
                        f_df.ix[which_forecast, which_force] = str((vdt1-vdt0).days) + 'd'
                    except:
                        pass
            else:
                f_df.ix[which_forecast, which_force] = 'no'
        except:
            # assume the directory is missing
            f_df.ix[which_forecast, which_force] = '--'

# what forecasts have been run successfully
r_dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
try:
    for item in os.listdir(r_dir0):
        if item[0] == 'f' and len(item) == 11:
            f_string = item
            fl = os.listdir(r_dir0 + f_string)

            if 'liveocean.in' in fl:
                f_df.ix[f_string, 'dot_in'] = 'YES'
            if 'low_passed.nc' in fl:
                f_df.ix[f_string, 'lp'] = 'YES'
            flh = [x for x in fl if 'ocean_his' in x]
            f_df.ix[f_string, 'his'] = str(int(flh[-1][10:14]))
except:
    pass
    
# what has been pushed to Azure (just the last num_days)
from azure.storage.blob import BlobService
azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
account = azu_dict['account']
key = azu_dict['key']
blob_service = BlobService(account_name=account, account_key=key)
for f_string in f_df.index[-args.num_days:]:
    ff_string = f_string.replace('.','')
    containername = ff_string
    try:
        blob_service.create_container(containername)
        blob_service.set_container_acl(containername, x_ms_blob_public_access='container')
        blobs = blob_service.list_blobs(containername)
        his_list = []
        for blob in blobs:
            his_list.append(blob.name)
        f_df.ix[f_string, 'azu'] = str(int(his_list[-1][-6:-3]))
    except:
        pass

# mark missing things
f_df[f_df.isnull()] = '--'

# make sure that the dates are in order
f_df = f_df.sort_index()

# print to the screen
#
# This prints beginning and end rows
#pd.set_option('display.max_rows', args.num_days)
#
# This prints just end rows
print(f_df[-args.num_days:])
