"""
Code to add dye to a portion of the domain of an existing history file.
"""

# imports
import numpy as np
import pickle
import netCDF4 as nc
import shutil
from time import time
import matplotlib.pyplot as plt

import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zrfun
import zfun

sys.path.append(os.path.abspath('../plotting'))
import pfun

# **** USER CHOICES *************
ex_name = 'lo8b'
ex_name_out = 'lo8dye'
f_string = 'f2019.07.03'
his_fn = 'ocean_his_0025.nc'
# *******************************

# Get Ldir
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['ex_name'] = ex_name
Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']

indir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
indir = indir0 + f_string + '/'
in_fn = indir + his_fn

outdir0 = indir0.replace(Ldir['ex_name'], ex_name_out)
outdir = outdir0 + f_string + '/'
out_fn = outdir + his_fn
Lfun.make_dir(outdir0)
Lfun.make_dir(outdir, clean=True)

# load dict of TEF segment indices
voldir = Ldir['LOo'] + 'tef2/' + 'volumes_' + Ldir['gridname'] + '/'
ji_dict = pickle.load(open(voldir + 'ji_dict.p', 'rb'))

tt0 = time()
# copy the file
shutil.copyfile(in_fn, out_fn)
print('Time to copy file = %0.1f sec' % (time()-tt0))

# get grid info
G = zrfun.get_basic_info(out_fn, only_G=True)

# create the dye field
ds = nc.Dataset(out_fn)
dye = 0 * ds['salt'][:]
ds.close()

seg_list = ['H'+str(item) for item in range(1,9)]

j_dict = {}; i_dict = {}
for seg_name in seg_list:
    jj = []; ii = []
    ji_list_full = ji_dict[seg_name]
    for ji in ji_list_full:
        jj.append(ji[0])
        ii.append(ji[1])
    jjj = np.array(jj)
    iii = np.array(ii)
    j_dict[seg_name] = jjj
    i_dict[seg_name] = iii
    
for seg_name in seg_list:    
    jjj = j_dict[seg_name]
    iii = i_dict[seg_name]
    dye[:,:,jjj,iii] = 1.0

foo = nc.Dataset(out_fn, 'a')
vv = foo.createVariable('dye_01', float, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
vv.long_name = 'passive dye'
vv.units = 'kg m-3'
vv[:] = dye
foo.close()

if Ldir['env'] == 'pm_mac':
    # PLOTTING
    plt.close('all')
    ds = nc.Dataset(out_fn)
    d1 = ds['dye_01'][0,-1,:,:]

    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111)

    cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], d1[1:-1,1:-1])
    pfun.dar(ax)
    pfun.add_coast(ax)
    fig.colorbar(cs)

    plt.show()
