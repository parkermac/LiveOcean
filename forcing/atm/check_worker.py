"""
Check the output of make_forcing_worker.py.
"""
# specify which forcing this code is for
frc = 'atm'

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(alp)
import zfun

#Info = Lfun.csv_to_dict(Ldir['LOo'] + 'current_Info/'
#    + frc + '/Info_for_main.csv')
    
# override date
Info = dict()
Info['date_string'] = '2013.03.19'
Info['f_string'] = 'f' + Info['date_string']
    
# define the output location

odir = Ldir['LOo'] + Ldir['gtag'] +'/' + Info['f_string'] +'/' + frc + '/'
odl = os.listdir(odir)
odl.sort()
import netCDF4 as nc    
for nn in odl:
    if '.nc' in nn:
        print('=== ' + nn + ' ===')

        ds = nc.Dataset(odir + nn)
        for vv in ds.variables:
            print(ds.variables[vv])
            
var_dict = {'lwrad_down':'lrf_time',
        'Pair':'pair_time',
        'Qair':'qair_time',
        'rain':'rain_time',
        'swrad':'srf_time',
        'Tair':'tair_time',
        'Uwind':'wind_time',
        'Vwind':'wind_time'}
        
# PLOTTING
import matplotlib.pyplot as plt
plt.close()

NR = 3
NC = 3
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

rlist = [0,0,0,1,1,1,2,2,2]
clist = [0,1,2,0,1,2,0,1,2]
jj = 0
for vn in var_dict.keys():
    tn = var_dict[vn]
    ds = nc.Dataset(odir + vn + '.nc')
    t = ds.variables[tn][:]
    td = (t - t[0])/86400.
    v = ds.variables[vn][:, 10, 10]
    
    nnr = rlist[jj]
    nnc = clist[jj]
    ax = axes[nnr, nnc]
    
    ax.plot(td, v)
    ax.set_title(vn)
    
    jj += 1

plt.show()