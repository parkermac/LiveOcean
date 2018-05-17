"""
Code to plot a single day of hycom fields processed
by ocn1/make_forcing_main.py
"""

import os
import sys
import pickle
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import pfun

pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cas3', tag='v1')
import zfun
import zrfun

pth = os.path.abspath('../forcing/hycom1')
if pth not in sys.path:
    sys.path.append(pth)

pth = os.path.abspath('../forcing/ocn1')
if pth not in sys.path:
    sys.path.append(pth)
import Ofun

date_string = '2017.01.01'
frc = 'ocn1'

Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
Ldir['LOogf'] = (Ldir['LOog'] + 'f' + date_string + '/')
Ldir['LOogf_f'] = (Ldir['LOogf'] + frc + '/')
Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')

#========================================================================

this_mo = 1 # target month to get data from
this_iz = -1 # index of HyCOM depth level to use (-1 = top)

#========================================================================

# +++ load ecology CTD cast data +++

dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
year = 2017
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# limit the stations used, if desired
sta_list = [s for s in sta_df.index]# if ('WPA' not in s) and ('GYS' not in s)]
sta_df = sta_df.loc[sta_list,['Max_Depth', 'Latitude', 'Longitude']]

#========================================================================

# get the Ecology cast data

Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')

Cast_dict = dict()

for station in sta_list:
    #print(' - reading: ' + station)           
    
    casts = Casts[Casts['Station'] == station]   
    
    casts = casts.set_index('Date')    
    # casts['Z'] = -casts['Depth'] # and make a Z column
    casts = casts.loc[:,['Salinity', 'Temperature','Z']] # keep only selected columns
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)

    # get the CTD cast data for this station, in the nearest month
    cdv = castdates.month.values
    if len(cdv) > 0:
        imo = zfun.find_nearest_ind(cdv, this_mo)
        new_mo = cdv[imo]
        cast = casts[casts.index.month==new_mo]
    
        Cast = cast.set_index('Z') # reorganize so that the index is Z
        Cast = Cast.dropna() # clean up
    
        # store cast in a dict
        Cast_dict[station] = Cast
    
        # save the month
        sta_df.loc[station,'Month'] = new_mo

    else:
        print(station + ': no data')
        # mask = sta_df.index != station
        # sta_df = sta_df.loc[mask,:]
    
sta_df = sta_df.dropna()
sta_list = list(sta_df.index)
#========================================================================

# +++ get the processed HYCOM data for the chosen day +++

in_dir = Ldir['LOogf_fd']
# coordinates
xyz = pickle.load(open(in_dir + 'coord_dict.p', 'rb'))
x = xyz['lon']
y = xyz['lat']
z = xyz['z'] # packed bottom to top
# fields
fh = pickle.load(open(in_dir + 'fh' + date_string +'.p', 'rb'))
xfh = pickle.load(open(in_dir + 'xfh' + date_string +'.p', 'rb'))

# make arrays of distance from the grid center
X, Y = np.meshgrid(x,y)
XX, YY = zfun.ll2xy(X, Y, x.mean(), y.mean())

# pull out fields for a single z level
t = fh['t3d'][this_iz,:,:]
xt = xfh['t3d'][this_iz,:,:]
s = fh['s3d'][this_iz,:,:]
xs = xfh['s3d'][this_iz,:,:]

if not (t.mask == s.mask).all():
    print('mask error!')

# initialize the filled fields
tnew = t.copy()
snew = s.copy()

# vectors or 2-column arrays of the good points
# to feed to cKDTree
xyorig = np.array((XX[~t.mask],YY[~t.mask])).T
torig = t[~t.mask]
sorig = s[~t.mask]

#========================================================================

# append CTD values so the good points from HYCOM
x_sta = sta_df.loc[:,'Longitude'].values
y_sta = sta_df.loc[:,'Latitude'].values
xx_sta, yy_sta = zfun.ll2xy(x_sta, y_sta, x.mean(), y.mean())
xy_sta = np.stack((xx_sta,yy_sta), axis=1)

xyorig = np.concatenate((xyorig, xy_sta))

temp_list = []
salt_list = []
for station in Cast_dict.keys():
    
    Cast = Cast_dict[station]
    cz = Cast.index.values
    izc = zfun.find_nearest_ind(cz, z[this_iz])
    this_temp = Cast.iloc[izc]['Temperature']
    this_salt = Cast.iloc[izc]['Salinity']
    temp_list.append(this_temp)
    salt_list.append(this_salt)
    
    sta_df.loc[station,'temp'] = this_temp
    sta_df.loc[station,'salt'] = this_salt


torig = np.concatenate((torig, np.array(temp_list,ndmin=1)))
sorig = np.concatenate((sorig, np.array(salt_list,ndmin=1)))

# array of the missing points that we want to fill
xynew = np.array((XX[t.mask],YY[t.mask])).T

# array of indices for all the missing points
a = cKDTree(xyorig).query(xynew)
aa = a[1]

# use those indices to fill in from the good data
tnew[t.mask] = torig[aa]
snew[t.mask] = sorig[aa]

# plotting
plt.close('all')
fig, axes = plt.subplots(1,2, squeeze=False, figsize=(12,7))

v0 = 4
v1 = 9

ax = axes[0,0]
cs = ax.pcolormesh(x,y,t, cmap='rainbow', vmin=v0, vmax=v1)
#fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('HYCOM original')

ax = axes[0,1]
cs = ax.pcolormesh(x,y,tnew, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('extrap. with CTD')
# add station locations
for station in sta_list:
    sta_x = sta_df.loc[station,'Longitude']
    sta_y = sta_df.loc[station,'Latitude']
    sta_temp = sta_df.loc[station,'temp']
    ax.text(sta_x, sta_y, '%0.1f' % (sta_temp),
    horizontalalignment='center',
    verticalalignment='center')

plt.show()


