"""
TEF functions.
"""
import pandas as pd
import netCDF4 as nc
import numpy as np

# path to alpha provided by driver
import zfun
import zrfun

def get_sect_df():
    # section definitions
    # * x and y are latitude and longitude and we require sections to be NS or EW so
    # either x0=x1 or y0=y1
    # * landward is the sign to multipy the transport by to be landward (1 or -1)
    sect_df = pd.DataFrame(columns=['x0', 'x1', 'y0', 'y1', 'landward', 'long_name'])
    
    sect_df.loc['JdFmouth',:] = [-124.6, -124.6, 48.35, 48.6, 1, 'Juan de Fuca Mouth']
    sect_df.loc['SoGnorth',:] = [-125.4, -124.6, 50, 50, -1, 'Strait of Georgia North']
    sect_df.loc['AInorth',:] = [-122.78, -122.58, 48.13, 48.13, -1, 'Admiralty Inlet North']
    
    return sect_df
    
def get_inds(x0, x1, y0, y1, G):
    
    # determine the direction of the section
    if (x0==x1) and (y0!=y1):
        sdir = 'NS'
    elif (x0!=x1) and (y0==y1):
        sdir = 'EW'
    else:
        print('Input points do not form a proper section')
        sdir='bad'
        sys.exit()
    
    # we assume a plaid grid, as usual
    if sdir == 'NS':
        lon = G['lon_u'][0,:].squeeze()
        lat = G['lat_u'][:,0].squeeze()
    elif sdir == 'EW':
        lon = G['lon_v'][0,:].squeeze()
        lat = G['lat_v'][:,0].squeeze()
        
    # we get all 4 i's or j's but only 3 are used
    i0, i1, fr = zfun.get_interpolant(np.array([x0]), lon, extrap_nan=True)
    if np.isnan(fr):
        print('Bad x point')
        sys.exit()
    else:
        ii0 = int(i0)
    i0, i1, fr = zfun.get_interpolant(np.array([x1]), lon, extrap_nan=True)
    if np.isnan(fr):
        print('Bad x point')
        sys.exit()
    else:
        ii1 = int(i1)
    j0, j1, fr = zfun.get_interpolant(np.array([y0]), lat, extrap_nan=True)
    if np.isnan(fr):
        print('Bad y0 point')
        sys.exit()
    else:
        jj0 = int(j0)
    j0, j1, fr = zfun.get_interpolant(np.array([y1]), lat, extrap_nan=True)
    if np.isnan(fr):
        print('Bad y1 point')
        sys.exit()
    else:
        jj1 = int(j1)

    # get mask and trim indices
    # Note: the mask in G is True on water points
    if sdir == 'NS':
        mask = G['mask_u'][jj0:jj1+1, ii0]
        # Note: argmax finds the index of the first True in this case
        igood0 = np.argmax(mask)
        igood1 = np.argmax(mask[::-1])
        # keep one mask point on each end, just to be sure we have a closed section
        Mask = mask[igood0-1:-igood1+1]
        # and change the indices to match.  These will be the indices
        # of the start and end points.
        jj0 = jj0 + igood0 - 1
        jj1 = jj1 - igood1 + 1
        print('  sdir=%2s: jj0=%4d, jj1=%4d, ii0=%4d' % (sdir, jj0, jj1, ii0))
        Lat = lat[jj0:jj1+1]
        Lon = lon[ii0] * np.ones_like(Mask)
    elif sdir == 'EW':
        mask = G['mask_v'][jj0, ii0:ii1+1]
        igood0 = np.argmax(mask)
        igood1 = np.argmax(mask[::-1])
        Mask = mask[igood0-1:-igood1+1]
        ii0 = ii0 + igood0 - 1
        ii1 = ii1 - igood1 + 1
        print('  sdir=%2s: jj0=%4d, ii0=%4d, ii1=%4d' % (sdir, jj0, ii0, ii1))
        Lon = lon[ii0:ii1+1]
        Lat = lat[jj0] * np.ones_like(Mask)
        
    return ii0, ii1, jj0, jj1, sdir, Lon, Lat, Mask
    
def start_netcdf(fn, out_fn, NT, NX, NZ, Lon, Lat, Ldir):
    # generating some lists
    vn_list = []
    ds = nc.Dataset(fn)
    if False:
        # all 3D variables on the s_rho grid
        for vv in ds.variables:
            vdim = ds.variables[vv].dimensions
            if ( ('ocean_time' in vdim) and ('s_rho' in vdim) ):
                vn_list.append(vv)
    else:
        # override
        vn_list.append('salt')
        vn_list.append('temp')
        # vn_list.append('NO3')
        # vn_list.append('oxygen')
    # and some dicts of long names and units
    long_name_dict = dict()
    units_dict = dict()
    for vn in vn_list + ['ocean_time']:
        try:
            long_name_dict[vn] = ds.variables[vn].long_name
        except:
            long_name_dict[vn] = ''
        try:
            units_dict[vn] = ds.variables[vn].units
        except:
            units_dict[vn] = ''
    ds.close()
    # add custom dict fields
    long_name_dict['q'] = 'transport'
    units_dict['q'] = 'm3 s-1'
    long_name_dict['lon'] = 'longitude'
    units_dict['lon'] = 'degrees'
    long_name_dict['lat'] = 'latitude'
    units_dict['lat'] = 'degrees'

    # initialize netcdf output file
    foo = nc.Dataset(out_fn, 'w')
    foo.createDimension('xi_sect', NX)
    foo.createDimension('s_rho', NZ)
    foo.createDimension('ocean_time', NT)
    foo.createDimension('sdir_str', 2)
    for vv in ['ocean_time']:
        v_var = foo.createVariable(vv, float, ('ocean_time',))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in vn_list + ['q']:
        v_var = foo.createVariable(vv, float, ('ocean_time', 's_rho', 'xi_sect'))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in ['lon', 'lat']:
        v_var = foo.createVariable(vv, float, ('xi_sect'))
        v_var.long_name = long_name_dict[vv]
        v_var.units = units_dict[vv]
    for vv in ['zeta']:
        v_var = foo.createVariable(vv, float, ('ocean_time', 'xi_sect'))
        v_var.long_name = 'Free Surface Height'
        v_var.units = 'm'

    # add static variables
    foo['lon'][:] = Lon
    foo['lat'][:] = Lat

    # add global attributes
    foo.gtagex = Ldir['gtagex']
    foo.date_string0 = Ldir['date_string0']
    foo.date_string1 = Ldir['date_string1']

    foo.close()
    
    return vn_list
    
def add_fields(ds, count, vn_list, G, S, sinfo):
    
    ii0, ii1, jj0, jj1, sdir, landward, NT, NX, NZ, out_fn = sinfo
    
    foo = nc.Dataset(out_fn, 'a')
    
    # get depth and dz
    if sdir=='NS':
        h = ds['h'][jj0:jj1+1,ii0:ii1+1].squeeze()
        zeta = ds['zeta'][0,jj0:jj1+1,ii0:ii1+1].squeeze()
        z = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(z, axis=0)
        DZ = dz.mean(axis=2)
        dy = G['DY'][jj0:jj1+1,ii0:ii1+1].squeeze()
        DY = dy.mean(axis=1)
        zeta = zeta.mean(axis=1)
    elif sdir=='EW':
        h = ds['h'][jj0:jj1+1,ii0:ii1+1].squeeze()
        zeta = ds['zeta'][0,jj0:jj1+1,ii0:ii1+1].squeeze()
        z = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(z, axis=0)
        DZ = dz.mean(axis=1)
        dy = G['DY'][jj0:jj1+1,ii0:ii1+1].squeeze()
        DY = dy.mean(axis=0)
        zeta = zeta.mean(axis=0)
    # and then create the array of cell areas on the section
    DA = DY.reshape((1, NX)) * DZ
    # then velocity and hence transport
    if sdir=='NS':
        vel = ds['u'][0, :, jj0:jj1+1, ii0].squeeze()
    elif sdir=='EW':
        vel = ds['v'][0, :, jj0, ii0:ii1+1].squeeze()
    q = vel * DA * landward
    
    foo['q'][count, :, :] = q
    foo['zeta'][count, :] = zeta
    foo['ocean_time'][count] = ds['ocean_time'][0]
    
    # save the tracer fields averaged onto this section
    for vn in vn_list:
        if sdir=='NS':
            vvv = (ds[vn][0,:,jj0:jj1+1,ii0].squeeze()
                + ds[vn][0,:,jj0:jj1+1,ii1].squeeze())/2
        elif sdir=='EW':
            vvv = (ds[vn][0,:,jj0,ii0:ii1+1].squeeze()
                + ds[vn][0,:,jj1,ii0:ii1+1].squeeze())/2
        foo[vn][count,:,:] = vvv
        
    foo.close()
    
