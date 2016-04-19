"""
This is the main program for making the LOW PASS forcing file.

Also this is the first attempt to use the new forcing functions.

Performance: 44 sec per day (mac) or ~5 hours per year of days.
    But, it took 70 sec on fjord, and 12 hours for a year
    (no bio variables).
"""

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import forcing_functions as ffun
Ldir, Lfun = ffun.intro()

# hack 4/11/2016 to handle the location of the new lobio1 output
#if ((Ldir['roms'] == '/pmr1/parker/LiveOcean_roms/')
#    & (Ldir['gtagex'] == 'cascadia1_base_lobio1')):
#    Ldir['roms'] = '/pmr2/darr/LiveOcean_roms/'
#print(Ldir['roms'])
# ****************** CASE-SPECIFIC CODE *****************
import zfun

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
    if Ldir['run_type'] == 'backfill':
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
        outfile = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
            '/f' + Ldir['date_string'] + '/low_passed.nc')
    elif Ldir['run_type'] == 'forecast':
        # use the middle day of the last forecast (= yesterday)
        # and today and tomorrow from today's forecast
        date_string = Ldir['date_string']
        from datetime import datetime, timedelta
        dt_now = datetime.strptime(Ldir['date_string'], '%Y.%m.%d')
        dt_yesterday = dt_now - timedelta(1)
        dt_list = [dt_yesterday, dt_now]
        for dt in dt_list:
            if dt == dt_yesterday:
                date_string = dt.strftime(format='%Y.%m.%d')
                indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                    '/f' + date_string + '/')
                for ii in range(2,26): # use range(2,26) to use Godin 71 hour filter
                    hnum = ('0000' + str(ii))[-4:]
                    flist.append(indir + 'ocean_his_' + hnum + '.nc')
            elif dt == dt_now:
                date_string = dt.strftime(format='%Y.%m.%d')
                indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                    '/f' + date_string + '/')
                for ii in range(2,49): # use range(2,49) to use Godin 71 hour filter
                    hnum = ('0000' + str(ii))[-4:]
                    flist.append(indir + 'ocean_his_' + hnum + '.nc')
        # make output name (full path)
        outfile = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
            '/f' + Ldir['date_string'] + '/low_passed.nc')
            
        # Old code that made tomorrow's low pass.  This may come in handy
        # when we start nesting.
        #indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
        #    '/f' + Ldir['date_string'] + '/')
        #for ii in range(2,73): # for Godin 71 hour filter
        #    hnum = ('0000' + str(ii))[-4:]
        #    flist.append(indir + 'ocean_his_' + hnum + '.nc')
        ## make output name (full path)
        #outfile = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
        #    '/f' + Ldir['date_string'] + '/low_passed_tomorrow.nc')
        
    # RUN THE FUNCTION
    roms_low_pass(flist, outfile, zfun)
    
    result_dict['result'] = 'success'    
except:
    result_dict['result'] = 'fail'
    
# ************** END CASE-SPECIFIC CODE *****************

ffun.finale(result_dict, Ldir, Lfun)

