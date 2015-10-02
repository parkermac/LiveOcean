"""
This is the main program for making the LOW PASS forcing file.

Also this is the first attempt to use the new forcing functions.

Still have to edit driver_forcing.sh.
"""

import os; import sys; fpth = os.path.abspath('../')
if fpth not in sys.path: sys.path.append(fpth)
import forcing_functions as ffun; reload(ffun)
Ldir, Lfun = ffun.intro()
# I think we can get zfun because ffun.intro() added alpha to the path.
import zfun; reload(zfun)

# ****************** CASE-SPECIFIC CODE *****************
result_dict = dict()   
try:
    
    # define the filtering function
    def roms_low_pass(flist, outfile, zfun):
        # create the filter
        nf = len(flist)
        if nf == 71:
            print(' - Using Godin filter')
            filt0 = zfun.godin_shape()
        else:
            print(' - Using Hanning filter for list length = ' + str(nf))
            filt0 = zfun.hanning_shape(nf)
        # create the output file
        import shutil
        shutil.copyfile(flist[0],outfile)
        # create the Datasets
        import netCDF4 as nc               
        ds = nc.MFDataset(flist)
        dsout = nc.Dataset(outfile,'a')
        # loop over all variables that have time axes
        for vn in ds.variables:
            if 'ocean_time' in ds.variables[vn].dimensions:
                #print(vn + ' ' + str(ds.variables[vn].shape)) # debugging
                ndim = len(ds.variables[vn].shape)
                filt_shape = (nf,)
                for ii in range(ndim-1):
                    filt_shape = filt_shape + (1,)
                v = ds.variables[vn][:]
                filt = filt0.reshape(filt_shape)
                vf = (filt*v).sum(axis=0)
                dsout.variables[vn][:] = vf.reshape(dsout.variables[vn].shape)
        ds.close()
        dsout.close()
        
    # make input list (full paths)
    flist = []
    # create the list of history files
    date_string = Ldir['date_string']
    from datetime import datetime, timedelta
    dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
    dt_tomorrow = dt_now + timedelta(1)
    dt_yesterday = dt_now - timedelta(1)
    dt_list = [dt_yesterday, dt_now, dt_tomorrow]
    for dt in dt_list:
        date_string = dt.strftime(format='%Y.%m.%d')
        indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
            '/f' + date_string + '/')
        for ii in range(2,26): # use range(2,26) to use Godin 71 hour filter
            hnum = ('0000' + str(ii))[-4:]
            flist.append(indir + 'ocean_his_' + hnum + '.nc')
    # remove the last item
    flist.pop() # cute
    # make output name (full path)
    outfile = indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
        '/f' + Ldir['date_string'] + '/low_passed.nc')
     
    # RUN THE FUNCTION
    roms_low_pass(flist, outfile, zfun)
    
    result_dict['result'] = 'success'    
except:
    result_dict['result'] = 'fail'
    
# ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)

