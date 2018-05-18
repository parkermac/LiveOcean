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

# User input: eventually any of these might be values in a list
# that is looped over

year = 2017 # target year
this_mo = 1 # target month to get data from
this_iz = 0 # index of HyCOM depth level to use (0=bottom to 39=top, or -1=top)

#========================================================================

def checknan(fld):
    if np.isnan(fld).sum() > 0:
        print('WARNING: nans in data field')    

#========================================================================

# +++ load ecology CTD cast data +++

dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# limit the stations used, if desired
sta_list = [s for s in sta_df.index]# if ('WPA' not in s) and ('GYS' not in s)]
# keep only certain columns
sta_df = sta_df.loc[sta_list,['Max_Depth', 'Latitude', 'Longitude']]

#========================================================================

# +++ load the Ecology cast data +++

Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')

# start a dict to store one cast per station (it is has data in the year)
Cast_dict = dict()

for station in sta_list:
    casts = Casts[Casts['Station'] == station]
    casts = casts.set_index('Date')
    casts = casts.loc[:,['Salinity', 'Temperature','Z']] # keep only selected columns
    # identify a single cast by its date
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    
    # get the CTD cast data for this station, in the nearest month
    cdv = castdates.month.values # all the months with casts
    if len(cdv) > 0:
        # get the cast closest to the selected month
        imo = zfun.find_nearest_ind(cdv, this_mo)
        new_mo = cdv[imo]
        cast = casts[casts.index==castdates[imo]]
        Cast = cast.set_index('Z') # reorganize so that the index is Z
        Cast = Cast.dropna() # clean up
        # store cast in a dict
        Cast_dict[station] = Cast
        # save the month, just so we know
        sta_df.loc[station,'Month'] = new_mo
    else:
        print(station + ': no data')
        
# at this point Cast_dict.keys() is the "official" list of stations
# to loop over

#========================================================================

# +++ get the processed HYCOM data for the chosen day +++

in_dir = Ldir['LOogf_fd']
# coordinates
xyz = pickle.load(open(in_dir + 'coord_dict.p', 'rb'))
x = xyz['lon']
y = xyz['lat']
z = xyz['z'] # packed bottom to top
zz = z[this_iz] # HYCOM target depth

# fields
fh = pickle.load(open(in_dir + 'fh' + date_string +'.p', 'rb'))

# make arrays of distance from the grid center
X, Y = np.meshgrid(x,y)
XX, YY = zfun.ll2xy(X, Y, x.mean(), y.mean())

# pull out fields for a single z level
t = fh['t3d'][this_iz,:,:]
s = fh['s3d'][this_iz,:,:]

# fill values to use for levels with all missing data
t0 = fh['t3d'].min()
s0 = fh['s3d'].max()

# and print infor about their locations
for vn in ['t3d', 's3d']:
    [k,m,l] = np.unravel_index(fh[vn].argmin(), fh[vn].shape)
    print('%s min at [%d, %d, %d]' % (vn, k, l, m))
    [k,m,l] = np.unravel_index(fh[vn].argmax(), fh[vn].shape)
    print('%s max at [%d, %d, %d]' % (vn, k, l, m))
                    
# check that the t.mask is the same as the s.mask
if not (t.mask == s.mask).all():
    print('mask error!')

# vectors or 1- or 2-column arrays of the good points to feed to cKDTree
xyorig = np.array((XX[~t.mask],YY[~t.mask])).T
torig = t[~t.mask]
sorig = s[~t.mask]

#========================================================================

# +++ append good points from CTD data to our arrays +++

goodcount = 0

for station in Cast_dict.keys():
    
    Cast = Cast_dict[station]
    cz = Cast.index.values
    izc = zfun.find_nearest_ind(cz, zz)
    
    # only take data from this cast if its bottom depth is at or above
    # the chosen hycom level
    czbot = -sta_df.loc[station,'Max_Depth']
    if czbot <= zz:
        # becasue we used find_nearest above we should always
        # get data in the steps below
        this_temp = Cast.iloc[izc]['Temperature']
        this_salt = Cast.iloc[izc]['Salinity']
        # and store in sta_df (to align with lat, lon)
        sta_df.loc[station,'temp'] = this_temp
        sta_df.loc[station,'salt'] = this_salt
        goodcount += 1
    else:
        pass
        
if goodcount >= 1:
    
    # drop stations that don't have T and s values at this depth
    sta_df = sta_df.dropna()
    # and for later convenience make a new list of stations
    sta_list = list(sta_df.index)
    
    # if we got any good points then append them
    print('goodcount = ' + str(goodcount) + ' =?')
    print('len(sta_df) = ' + str(len(sta_df)))
    
    # append CTD values to the good points from HYCOM
    x_sta = sta_df['Longitude'].values
    y_sta = sta_df['Latitude'].values
    xx_sta, yy_sta = zfun.ll2xy(x_sta, y_sta, x.mean(), y.mean())
    xy_sta = np.stack((xx_sta,yy_sta), axis=1)
    xyorig = np.concatenate((xyorig, xy_sta))
    
    temp_arr = sta_df['temp'].values
    salt_arr = sta_df['salt'].values
    torig = np.concatenate((torig, np.array(temp_arr,ndmin=1)))
    sorig = np.concatenate((sorig, np.array(salt_arr,ndmin=1)))
    
else:
    print('No points added')
    
#========================================================================

# +++ do the extrapolation +++

def extrap_nearest_to_nan_CTD(XX,YY,fld,xyorig=[],fldorig=[],fld0=0):
    if np.ma.is_masked(fld):
        if fld.all() is np.ma.masked:
            print('  filling with ' + str(fld0))
            fldf = fld0 * np.ones(fld.data.shape)
            fldd = fldf.data
            checknan(fldd)
            return fldd
        else:
            fldf = fld.copy()
            # array of the missing points that we want to fill
            xynew = np.array((XX[fld.mask],YY[fld.mask])).T
            # array of indices for points nearest to the missing points
            a = cKDTree(xyorig).query(xynew)
            aa = a[1]

            # use those indices to fill in using the good data
            fldf[fld.mask] = fldorig[aa]
                
            fldd = fldf.data
            checknan(fldd)
            return fldd
    else:
        checknan(fld)
        return fld
        
tnew = extrap_nearest_to_nan_CTD(XX,YY,t,xyorig=xyorig,fldorig=torig,fld0=t0)
snew = extrap_nearest_to_nan_CTD(XX,YY,s,xyorig=xyorig,fldorig=sorig,fld0=s0)

# # initialize the fields to be filled
# tnew = t.copy()
# snew = s.copy()
# # array of the missing points that we want to fill
# xynew = np.array((XX[t.mask],YY[t.mask])).T
# # array of indices for points nearest to the missing points
# a = cKDTree(xyorig).query(xynew)
# aa = a[1]
# # use those indices to fill in using the good data
# tnew[t.mask] = torig[aa]
# snew[t.mask] = sorig[aa]

#========================================================================

# +++ plotting +++

plt.close('all')
fig, axes = plt.subplots(1,2, squeeze=False, figsize=(14,7))

# temperature color limts
v0 = 4
v1 = 9

# ax = axes[0,0]
# cs = ax.pcolormesh(x,y,t, cmap='rainbow', vmin=v0, vmax=v1)
# #fig.colorbar(cs, ax=ax)
# pfun.dar(ax)
# pfun.add_coast(ax)
# ax.set_xlim(x[0],x[-1])
# ax.set_ylim(y[0],y[-1])
# ax.set_title('HYCOM temp: n=' + str(this_iz) + ', z=' + str(zz))

ax = axes[0,0]
cs = ax.pcolormesh(x,y,tnew, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('temp w/CTD: n=' + str(this_iz) + ', z=' + str(zz))
# add station locations
if goodcount >= 1:
    for station in sta_list:
        sta_x = sta_df.loc[station,'Longitude']
        sta_y = sta_df.loc[station,'Latitude']
        sta_temp = sta_df.loc[station,'temp']
        ax.text(sta_x, sta_y, '%0.1f' % (sta_temp),
        horizontalalignment='center',
        verticalalignment='center')

# salt color limits
v0 = 20
v1 = 33

ax = axes[0,1]
cs = ax.pcolormesh(x,y,snew, cmap='rainbow', vmin=v0, vmax=v1)
fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('salt w/CTD')
# add station locations
if goodcount >= 1:
    for station in sta_list:
        sta_x = sta_df.loc[station,'Longitude']
        sta_y = sta_df.loc[station,'Latitude']
        sta_salt = sta_df.loc[station,'salt']
        ax.text(sta_x, sta_y, '%0.1f' % (sta_salt),
        horizontalalignment='center',
        verticalalignment='center')

plt.show()
    


