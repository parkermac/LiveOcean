"""
Functions for the surface code.
"""

# imports
import os
import netCDF4 as nc
import zfun

def create_ds(out_fn):
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    out_ds = nc.Dataset(out_fn, 'w')
    return out_ds
    
def add_dims(in_ds, out_ds):
    # lists of essential dimensions and variables
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
    # Copy dimensions
    for dname, the_dim in in_ds.dimensions.items():
        if dname in dlist:
            out_ds.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy variables
    for vn in vn_list2:
        varin = in_ds[vn]
        vv = out_ds.createVariable(vn, varin.dtype, varin.dimensions)
        vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        vv[:] = in_ds[vn][:]
        
def add_fields(in_ds, out_ds, vn_list2t, vn_list3t, slev=-1, suffix=''):
    # slev = 0 for bottom, -1 for top
    for vn in vn_list2t:
        varin = in_ds[vn]
        vv = out_ds.createVariable(vn, varin.dtype, varin.dimensions)
        vv.long_name = varin.long_name
        vv.units = varin.units
        try:
            vv.time = varin.time
        except AttributeError:
            # ocean_time has no time
            pass
        vv[:] = zfun.fillit(in_ds[vn][:])
    #
    for vn in vn_list3t:
        do_var = True
        try:
            varin = in_ds[vn]
        except IndexError:
            # designed so that it still works when we don't have a variable from this list
            # e.g. when there is no bio or carbon
            do_var = False
            print(' - Variable not found: ' + vn)
        if do_var==True:
            dd = tuple([d for d in varin.dimensions if d != 's_rho'])
            vv = out_ds.createVariable(vn + suffix, varin.dtype, dd)
            if vn == 'PH':
                vv.long_name = 'pH'
            elif vn == 'ARAG':
                vv.long_name = 'Aragonite Saturation State'
            else:
                vv.long_name = varin.long_name
            try:
                vv.units = varin.units
            except AttributeError:
                # salt has no units
                pass
            vv.time = varin.time
            vv[:] = zfun.fillit(in_ds[vn][:, slev, :, :].squeeze())
            
def add_carbon():
    """
    This code will eventually use the python-native co2sys code to add pH and Omega_arag to the output files.
    
    Notes from carbon_fun.py in LiveOcean/co2sys:
    
    variables we get: ['alkalinity', 'TIC', 'salt', 'temp','rho'] (into v_dict)
    
    # create depth, pressure
    depth = -z_rho[NZ, :, :].squeeze()
    pres = sw.pres(depth, lat)
    
    # assume potential temperature is close enough, or run this:
    # temp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, v_dict['pres'])
    
    # convert from umol/L to umol/kg using in situ dentity
    v_dict['alkalinity'] = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
    v_dict['TIC'] = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
    
    and here is what we do in matlab:
    A = CO2SYS(a.alkalinity(:), a.TIC(:), 1, 2, a.salt(:), a.temp(:), a.temp(:), a.pres(:), a.pres(:), 50, 2, 1, 10, 1);
    ph = A(:,18);
    om = A(:,31);
    PH = reshape(ph, size(a.salt));
    OM = reshape(om, size(a.salt));
    PH = real(PH);
    OM = real(OM);
    """
            
            