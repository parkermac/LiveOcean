"""
Plots the HYCOM data.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
import zfun; reload(zfun)
import matfun; reload(matfun)
import netCDF4 as nc

# get the coastline
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'
cmat = matfun.loadmat(fn_coast)

# PLOTTING
import matplotlib.pyplot as plt
cmap = plt.get_cmap(name='rainbow')
plt.close()

fig, axes = plt.subplots(nrows=3, ncols=5, figsize=(16, 12), squeeze=False)

vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']
tag = '_combined' # e.g. '91.0'
nc_dir = Ldir['data'] + 'hycom' + tag + '/'

#nc_dir = Ldir['LOo'] + Ldir['gtag'] + '/f2015.02.14/ocn/Data/'

cnum = 0
for vn in vn_list:
    
    fn = nc_dir + vn + '.nc'
    ds = nc.Dataset(fn)
    #zfun.ncd(ds)

    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    tmod = ds.variables['tmod'][:] # ['dt'] in the forecast version?
    
    #nt = int(len(tmod)/2)
    nt=-1
    # target position (-124.5, 47 = RN)
    Lon = lon[0,:]
    Lat = lat[:,0]
    xbool = Lon >= -124.5
    ybool = Lat >= 47
    ii = list(xbool).index(True)
    jj = list(ybool).index(True)
    if vn == 'ssh':
        fldz = 0      
        fld = ds.variables[vn][nt,:,:].squeeze()
        flde = ds.variables[vn + '_extrap'][nt,:,:].squeeze()
        fldf = ds.variables[vn + '_filt'][nt,:,:].squeeze()
        # get time series
        fs = ds.variables[vn][:, jj, ii].squeeze()
        fse = ds.variables[vn + '_extrap'][:, jj, ii].squeeze()
        fsf = ds.variables[vn + '_filt'][:, jj, ii].squeeze()
    else:
        kk = -1
        fldz = ds.variables['z'][kk]
        fld = ds.variables[vn][nt,kk,:,:].squeeze()
        flde = ds.variables[vn + '_extrap'][nt,kk,:,:].squeeze()
        fldf = ds.variables[vn + '_filt'][nt,kk,:,:].squeeze()
        # get time series
        fs = ds.variables[vn][:,kk, jj, ii].squeeze()
        fse = ds.variables[vn + '_extrap'][:, kk, jj, ii].squeeze()
        fsf = ds.variables[vn + '_filt'][:, kk, jj, ii].squeeze()
   
    ds.close()

    ax0 = axes[0, cnum]
    ax1 = axes[1, cnum]
    ax2 = axes[2, cnum]

    cs = ax0.pcolormesh(lon, lat, fld, cmap = cmap)
    ax0.plot(lon[jj, ii], lat[jj, ii],'*r')
    ax0.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    ax0.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    ax0.xaxis.set_ticklabels([])
    ax0.yaxis.set_ticklabels([])
    zfun.dar(ax0)
    ax0.set_title(vn)
    fig.colorbar(cs, ax=ax0)

    cs = ax1.pcolormesh(lon, lat, flde, cmap = cmap)
    ax1.plot(lon[jj, ii], lat[jj, ii],'*r')
    ax1.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    ax1.axis([lon.min(), lon.max(), lat.min(), lat.max()])
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    zfun.dar(ax1)
    fig.colorbar(cs, ax=ax1)
    
    from datetime import datetime
    tmod0 = Lfun.datetime_to_modtime(datetime(2013,1,1))
    tdays = (tmod - tmod0)/86400.
    ax2.plot(tdays, fs,'-*b', linewidth=3 )
    ax2.plot(tdays, fsf,'-or', linewidth=3 )
    ax2.set_title(vn + ' ' + str(fldz) + ' m')
    ax2.set_xlabel('Days from 1/1/2013')
    ax2.set_xlim(tdays[0], tdays[-1])
    # add year divisions
    aa = ax2.axis()
    ax2.plot([0,0],aa[2:],'-k')
    ax2.plot([365,365],aa[2:],'-k')
    ax2.plot([730,730],aa[2:],'-k')
        
    cnum += 1

plt.show()
plt.savefig('/Users/PM5/Desktop/hycom' + tag + '.png')


