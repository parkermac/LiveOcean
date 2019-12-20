"""
Plot results of a particle tracking experiment.
"""

# setup
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun
import zfun
import zrfun
sys.path.append(os.path.abspath('../plotting'))
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

Ldir = Lfun.Lstart()

# Choose an experiment to plot from.
indir0 = Ldir['LOo'] + 'tracks2/'
indir_list_raw = os.listdir(indir0)
indir_list = []
for d in indir_list_raw:
    if os.path.isdir(indir0 + d):
        indir_list.append(d)
indir_list.sort()
Npt = len(indir_list)#
print('\n%s\n' % '** Choose Experiment to plot **')
for npt in range(Npt):
    print(str(npt) + ': ' + indir_list[npt])
my_npt = input('-- Experiment number (return = 0) --')
if len(my_npt)==0:
    my_npt = 0
indir = indir_list[int(my_npt)] + '/'

# Choose a release from this experiment.
rel_list = [rel for rel in os.listdir(indir0 + indir) if 'release' in rel]
rel_list.sort()
Nrl = len(rel_list)
print('\n%s\n' % '** Choose Release file to plot **')
for nrl in range(Nrl):
    print(str(nrl) + ': ' + rel_list[nrl])
my_nrl = input('-- Release number (return = 0) -- ')
if len(my_nrl)==0:
    my_nrl = 0
rel = rel_list[int(my_nrl)]

# get Datasets
dsr = nc4.Dataset(indir0 + indir + rel)
dsg = nc4.Dataset(indir0 + indir + 'grid.nc')

# and a dict of run info
EI = Lfun.csv_to_dict(indir0 + indir + 'exp_info.csv')

# get grid info
G = zrfun.get_basic_info('/Users/pm7/Documents/LiveOcean_roms/output/' + 
    'cas6_v3_lo8b/f2019.07.04/ocean_his_0001.nc', only_G=True)
xvec = G['lon_rho'][0,:]
yvec = G['lat_rho'][:,0]

# The goal is to figure out how many of the particles do not lie in the
# region of the initial release.  To do this we will build lists of
# the (j,i) indices of the particles.

# First gather the ji list of the initial release
if EI['ic_name'] == 'hc1': # Hood Canal using TEF "segments"
    # NOTE: this one requires information about grid indices in certain regions of
    # the Salish Sea, and is very specific to the analysis in x_tef.
    # Not for general use.
    import pickle
    # select the indir
    import Lfun
    Ldir = Lfun.Lstart()
    indir0 = Ldir['LOo'] + 'tef/cas6_v3_lo8b_2017.01.01_2017.12.31/'
    indir = indir0 + 'flux/'
    # load data
    ji_dict = pickle.load(open(indir + 'ji_dict.p', 'rb'))
    ji_list = []
    seg_list = ['H3','H4','H5','H6','H7','H8']
    for seg_name in seg_list:
        ji_list = ji_list + ji_dict[seg_name]
    #ji_set = set(ji_list)

# now find what fraction of the initial condition that is in this list
# as a test of the method
lon = dsr['lon'][:]
lat = dsr['lat'][:]
ot_vec = dsr['ot'][:]
days = (ot_vec - ot_vec[0])/86400
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

NT, NP = lon.shape

day_list = []
frac_list = []
for tt in range(0,NT,100):
    x = lon[tt,:]; y = lat[tt,:]
    this_ji_list = []
    for p in range(NP):
        ii = (np.abs(xvec - x[p])).argmin()
        jj = (np.abs(yvec - y[p])).argmin()
        this_ji_list.append((jj,ii))
    #this_ji_set = set(this_ji_list)
    
    # find the number of particles still in the initial region
    incount = 0
    for ji in this_ji_list:
        if ji in ji_list:
            incount += 1
            
    day_list.append(days[tt])
    frac_list.append(incount/len(this_ji_list))
            
    #print('%0.1f days: Wrong Fraction in initial volume = %0.2f' % (days[tt],len(this_ji_set & ji_set)/len(this_ji_set)))
    print('%0.1f days: Right Fraction in initial volume = %0.2f' % (days[tt],incount/len(this_ji_list)))
    
plt.close('all')
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

day_vec = np.array(day_list)
frac_vec = np.array(frac_list)

ax.plot(day_vec, frac_vec, '-g', linewidth=3)
ax.grid(True)
ax.set_ylim(0,1)

plt.show()


