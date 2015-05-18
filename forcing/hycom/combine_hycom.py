"""
Code to combine NetCDF files from different HYCOM exnum's
"""
# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
import zfun; reload(zfun)
import hfun; reload(hfun)
import netCDF4 as nc
import shutil

exnum_list = ['90.9', '91.0', '91.1']
vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']

# make a clean output directory
nc_dir_new = Ldir['data'] + 'hycom_combined/'
Lfun.make_dir(nc_dir_new, clean=True)

# get a dictionary of time axes.
# and make a copy of the first exnum file into "nc_dir_new"

# repeat this loop for each variable
for vn in vn_list:
    print 'Working on variable ' + vn
    
    # first get a dictionary of time axes
    nex = 0 # counter for which experiment number
    tmod_dict = dict()   
    for exnum in exnum_list:   
        nc_dir = Ldir['data'] + 'hycom' + exnum + '/'
        fn = nc_dir + vn + '.nc'
        ds = nc.Dataset(fn)
        tmod_dict[nex] = ds.variables['tmod'][:]
        ds.close()
        # and if this is the first exnum, initialize the combined netcdf
        # file by copying        
        if nex == 0:
            fn_new = nc_dir_new + vn + '.nc'
            print '  Copying first exnum into ' + fn_new
            shutil.copyfile(fn,fn_new) 
        nex += 1      

    # now loop through the exnums again, appending data as needed
    for nex in range(1,len(exnum_list)):
        
        # get the new file ready to write to    
        fn_new = nc_dir_new + vn + '.nc'
        foo = nc.Dataset(fn_new, 'r+')
        current_tmod = foo.variables['tmod'][:]
        
        # find the index to start appending at       
        istart = list(current_tmod).index(tmod_dict[nex][0])
        
        exnum = exnum_list[nex]
        print ''
        print '  Appending data from exnum ' + exnum
        print '  istart = ' + str(istart)
        
        # account for the last four times in 90.9 being bad data
        if exnum_list[nex-1] == '90.9':
            istart = istart - 5
                         
        # get the file that we are going to write from
        nc_dir = Ldir['data'] + 'hycom' + exnum + '/'
        fn = nc_dir + vn + '.nc'
        ds = nc.Dataset(fn)
        # and get the data from it
        fld = ds.variables[vn][:]
        flde = ds.variables[vn + '_extrap'][:]
        fldt = ds.variables['tmod'][:]
        
        # append from this exnum to the new file
        foo.variables['tmod'][istart:] = fldt
        foo.variables[vn][istart:] = fld
        foo.variables[vn + '_extrap'][istart:] = flde
                
        foo.close()
        ds.close()
        
        

    

