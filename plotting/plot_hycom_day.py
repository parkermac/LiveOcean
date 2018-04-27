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

pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cas3', tag='v1')
import zfun
import zrfun

date_string = '2017.01.01'
frc = 'ocn1'

Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
Ldir['LOogf'] = (Ldir['LOog'] + 'f' + date_string + '/')
Ldir['LOogf_f'] = (Ldir['LOogf'] + frc + '/')
Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')


#========================================================================

this_mo = 1
this_zlev = -11 # index of HyCOM depth level to use

#========================================================================

dir0 = Ldir['parent'] + 'ptools_data/ecology/'

# load Ecology station location and depth info
sta_info_fn = dir0 + 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx'
sta_df = pd.read_excel(sta_info_fn)
sta_df = sta_df.set_index('Station')
# get locations in decimal degrees
for sta in sta_df.index:
    lat_str = sta_df.loc[sta, 'Lat_NAD83 (deg / dec_min)']
    lat_deg = float(lat_str.split()[0]) + float(lat_str.split()[1])/60
    sta_df.loc[sta,'Latitude'] = lat_deg
    #
    lon_str = sta_df.loc[sta, 'Long_NAD83 (deg / dec_min)']
    lon_deg = float(lon_str.split()[0]) + float(lon_str.split()[1])/60
    sta_df.loc[sta,'Longitude'] = -lon_deg    
sta_df.pop('Lat_NAD83 (deg / dec_min)')
sta_df.pop('Long_NAD83 (deg / dec_min)')

# get the Ecology cast data
fn = dir0 + 'ParkerMacCready2017CTDDataFeb2018.xlsx'
# read in the data (all stations, all casts)
all_casts = pd.read_excel(fn, sheet_name='2017Provisional_CTDResults',
                          parse_dates = ['Date'])
# data long names; we retain only these fields
data_long_names = ['Salinity', 'Temp','Z']
sta_list = [sta for sta in sta_df.index]
#

#========================================================================

# get the ROMS grid
fng = Ldir['grid'] + 'grid.nc'
dsg = nc.Dataset(fng)
xr = dsg['lon_rho'][:]
yr = dsg['lat_rho'][:]
mr = dsg['mask_rho'][:]

# get the ini field on the ROMS grid
fn = Ldir['LOogf_f'] + 'ocean_ini.nc'
ds = nc.Dataset(fn)
tr = ds['temp'][0, this_zlev, :, :].squeeze()
# mask on land
tr[mr==0] = np.nan

# get the processed hycom data that went into it
in_dir = Ldir['LOogf_fd']
# coordinates
xyz = pickle.load(open(in_dir + 'coord_dict.p', 'rb'))
x = xyz['lon']
y = xyz['lat']
z = xyz['z']
# fields
fh = pickle.load(open(in_dir + 'fh' + date_string +'.p', 'rb'))
xfh = pickle.load(open(in_dir + 'xfh' + date_string +'.p', 'rb'))
t = fh['t3d'][this_zlev,:,:]
xt = xfh['t3d'][this_zlev,:,:]
tt = t.copy()


# ***************************************************************
for station in sta_list:
    print(' - reading: ' + station)           
    casts = all_casts[all_casts['Station'] == station]   
    casts = casts.set_index('Date')    
    casts['Z'] = -casts['Depth'] # and make a Z column
    casts = casts[data_long_names] # keep only selected columns
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    cast_info = station + ': ' + sta_df.loc[station,'Descrip']
    Max_z = -float(sta_df.loc[station, 'Max_Depth'])

    # get the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        if imo == this_mo:
            cast = casts[casts.index==cdate]
            # drop repeat values (aiming for the first of a depth pair)
            zdf = np.diff(cast['Z'])
            zdf = np.concatenate((np.array([1.,]),zdf))
            mask = zdf != 0
            cast = cast[mask]
            cast = cast[:-5] # drop bottom values (sometimes bad)
            
            Cast = cast.set_index('Z')

            sta_x = sta_df.loc[station,'Longitude']
            sta_y = sta_df.loc[station,'Latitude']

            cz = Cast.index.values
            izc = zfun.find_nearest_ind(cz, z[this_zlev])
            sta_temp = Cast.iloc[izc]['Temp']

            # try adding a point and then extrapolating
            i0, i1, ifr = zfun.get_interpolant(np.array([sta_x]), x)
            j0, j1, jfr = zfun.get_interpolant(np.array([sta_y]), x)
            tt[j0,i0] = sta_temp # doing this sets the mask to False automatically

# ***************************************************************

X, Y = np.meshgrid(x,y)

XX, YY = zfun.ll2xy(X, Y, -124, 46)

TT = Ofun.extrap_nearest_to_masked(XX, YY, tt)

# plotting
plt.close('all')
fig, axes = plt.subplots(1,3, squeeze=False, figsize=(12,7))

v0 = 6
v1 = 12

ax = axes[0,0]
cs = ax.pcolormesh(x,y,t, cmap='rainbow', vmin=v0, vmax=v1)
#fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('HyCOM original')

ax = axes[0,1]
cs = ax.pcolormesh(x,y,xt, cmap='rainbow', vmin=v0, vmax=v1)
#fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('extrapolated')

ax = axes[0,2]
cs = ax.pcolormesh(x,y,TT, cmap='rainbow', vmin=v0, vmax=v1)
#fig.colorbar(cs, ax=ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_title('extrap. with CTD')

# ax = axes[0,3]
# cs = ax.pcolormesh(xr,yr,tr, cmap='rainbow', vmin=v0, vmax=v1)
# #fig.colorbar(cs, ax=ax)
# pfun.dar(ax)
# pfun.add_coast(ax)
# ax.set_xlim(x[0],x[-1])
# ax.set_ylim(y[0],y[-1])
# ax.set_title('Interp to ROMS')


plt.show()


