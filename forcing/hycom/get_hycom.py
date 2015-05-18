"""
This program downloads a regional subset of HYCOM output for a time period,
and puts it in a NetCDF file.

The only processing is that it makes the vertical packing bottom-to-top.

A good place to start for detailed info on HYCOM is here:
http://hycom.org/data/glbu0pt08/expt-91pt1

Workflow Notes:
    
* get_hycom.py gets the regridded HYCOM output in our region for a specified
time period, for a given "experiment,"
and packs them in separate NetCDF files: [ssh,t3d,s3d,u3d,v3d].nc,
e.g. in LiveOcean_data/hycom90.9/.

* extrapolate_hycom.py needs to be run for each "exnum" before
doing anything else, making new variables with the tag "_extrap".

* combine_hycom.py joins together the outputs of extrapolate_hycom.py
into [ssh,t3d,s3d,u3d,v3d].nc in LiveOcean_data/hycom_combined/.

* filter_hycom.py applies a 5-day Hanning filterto the results in hycom_combined,
creating variables with the tag "_filt".

* plot_hycom.py plots the results in e.g. hycom_combined:
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
from datetime import datetime, timedelta

# This section sets the datenum limits for the HYCOM extraction to get.

# current output of "hycom_availability.py"
#Working on 90.9
# time units = hours since 2000-01-01 00:00:00
# dt start = 2012-05-13 00:00:00
# dt end   = 2013-08-20 00:00:00
#
#Working on 91.0
# time units = hours since 2000-01-01 00:00:00
# dt start = 2013-08-17 00:00:00
# dt end   = 2014-04-08 00:00:00
#
#Working on 91.1
# time units = hours since 2000-01-01 00:00:00
# dt start = 2014-04-07 00:00:00
# dt end   = 2015-03-12 00:00:00

exnum = '91.1'

# This section gets the time axis of a given experiment, in datetime format,
# and then determines the indices, n0 and n1, into this list that correspond
# to our chosen time limits.

print ''
print 'Working on ' + exnum
fn = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_' + exnum
sys.stdout.flush()
ds = nc.Dataset(fn)
if False: # set to True to see what is in the files
    zfun.ncd(ds)
# get the time in a meaningful format
t_hycom = ds.variables['time'][:].squeeze()
tu = ds.variables['time'].units
print ' time units = ' + tu # should be 'hours since 2000-01-01 00:00:00'
t_origin = ds.variables['time'].time_origin
dt00 = datetime.strptime(t_origin, '%Y-%m-%d %H:%M:%S')
dt_list = [] # initialize a list
for tt in t_hycom:    
    dt_list.append(dt00 + timedelta(tt/24.))
print ' dt start = ' + datetime.strftime(dt_list[0], '%Y-%m-%d %H:%M:%S')
print ' dt end   = ' + datetime.strftime(dt_list[-1], '%Y-%m-%d %H:%M:%S')
sys.stdout.flush()
ds.close()

nc_dir = Ldir['data'] + 'hycom' + exnum + '/'

def get_dt_list_stored(nc_dir):
    dt_list_stored = []
    import netCDF4 as nc
    # find the times that are already in the file
    fn_stored = nc_dir + 'ssh.nc'
    try:
        ds = nc.Dataset(fn_stored)
        t_stored = ds.variables['tmod'][:].squeeze()
        dt_list_stored = []
        for tt in t_stored:
            dt_list_stored.append(Lfun.modtime_to_datetime(tt))
        ds.close()
    except:
        pass # possible that no files exist for this exnum
    return dt_list_stored
    
dt_list_stored = get_dt_list_stored(nc_dir)
if len(dt_list_stored) > 0:
    dt0 = dt_list_stored[-1] # the last stored time
else:
    dt0 = dt_list[0]
       
cc = 0 # a counter, for testing

Testing = False

while (dt0 < dt_list[-1]) and (cc < 2):
       
    try:
        nt0 = dt_list.index(dt0) # its index in the hycom database for exnum
        print dt0
        print ''
        print '*** Next Interval ***'
        print 'nt0 = ' + str(nt0)
    except:
        print 'error with dt0!'
        nt0 = -1

    ntmax = len(dt_list) - 1
    
    if Testing == True:
        ntimes_to_get = 3
    else:
        ntimes_to_get = 30
    # it actually gets one more than this,
    # but the first is a repeat
      
    nt1 = nt0 + ntimes_to_get
    if nt1 >= ntmax:
        nt1 = ntmax
    try:        
        dt1 = dt_list[nt1]
        print dt1
        print 'nt1 = ' + str(nt1)
    except:
        print 'error with dt1!'
        nt1 = -1
    
    # Here we get the HYCOM data from the server, and save all variables in
    # a dict.
    print 'Getting HYCOM data for ' + str(nt1-nt0+1) + ' times'
    sys.stdout.flush()
    print 'nt0 = ' + str(nt0) + ' nt1 = ' + str(nt1)
    import time
    tt0 = time.time()         
    out_dict = hfun.get_hycom_past(fn, nt0, nt1)
    print '  took ' + str(round(time.time() - tt0)) + ' seconds'
    sys.stdout.flush()
    
    # Now we add the time axis to that dict, both in datetime format (dt),
    # and model time format (t) = seconds since 1/1/1970.
    out_dict['dt'] = dt_list[nt0:nt1+1]
    dt = out_dict['dt']
    t = []
    for this_dt in dt:
        t.append(Lfun.datetime_to_modtime(this_dt))
    out_dict['t'] = t
    
    # Here we decide what index numbers this data will occupy in the output
    # NetCDF files (one per variable)
    NT0 = nt0
    NT1 = nt1
    
    # only clobbber the target directory if we are starting a new file
    if NT0 == 0:
        print 'making clean folder ' + nc_dir
        Lfun.make_dir(nc_dir, clean=True)
    
    print 'Writing HYCOM data to NetCDF'
    tt0 = time.time()   
    hfun.hycom_dict_to_netcdf(out_dict, nc_dir, NT0, NT1)
    print '  took ' + str(round(time.time() - tt0)) + ' seconds'
    
    dt_list_stored = get_dt_list_stored(nc_dir)
    dt0 = dt_list_stored[-1] # the last stored time
    
    sys.stdout.flush()
    
    if Testing == True:
        cc += 1

    
