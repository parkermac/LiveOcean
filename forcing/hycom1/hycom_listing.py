"""
Find information about what hycom archived files are available.

"""

# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

from importlib import reload
import hfun
reload(hfun)

import netCDF4 as nc
from datetime import datetime
import pickle

out_dir = Ldir['data'] + 'hycom1/'
Lfun.make_dir(out_dir, clean=False)

exnum_list = ['90.9', '91.0', '91.1', '91.2']

dt_dict = dict()
coords_dict = dict()

for exnum in exnum_list:
    print('\nWorking on ' + exnum)
    fn = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_' + exnum  
    ds = nc.Dataset(fn)
    
    coords = hfun.get_coordinates(ds)
    
    dt_list = hfun.get_dt_list(ds)
    print(' dt start = ' + datetime.strftime(dt_list[0], '%Y-%m-%d %H:%M:%S'))
    print(' dt end   = ' + datetime.strftime(dt_list[-1], '%Y-%m-%d %H:%M:%S'))
    sys.stdout.flush()
    
    # check for missing days
    td = (dt_list[-1] - dt_list[0]).days
    nd = len(dt_list)
    nmissing = td - nd
    print(' Have %d out of %d days (missing %d)' % (nd, td, nmissing))
    
    # print grid info
    print(' Indices %d:%d %d:%d' % (coords['i0'], coords['i1'], coords['j0'], coords['j1']))
    # in principle we should also check that the coordinates did not change between exnums,
    # but it appears they do not so I will ignore that for now.
    
    # save the datetimes and coordinates
    dt_dict[exnum] = dt_list
    coords_dict[exnum] = coords
    
    ds.close()

# Save the time vectors in a pickled dict.
dt_fn = out_dir + 'dt_dict.p'  
pickle.dump(dt_dict, open(dt_fn, 'wb'))
# load this with a call like:
# dt_dict = pickle.load(open(out_fn, 'rb'))

# also save coordinates
xyz_fn = out_dir + 'coords_dict.p'  
pickle.dump(coords_dict, open(xyz_fn, 'wb'))
