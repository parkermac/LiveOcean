"""
Time filtering of HYCOM NetCDF files.

Do this AFTER extraploating.
"""
# setup
import os; import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

ofp = os.path.abspath('../ocn')
if ofp not in sys.path:
    sys.path.append(ofp)
import Ofun; reload(Ofun)

vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']
tag = '_combined' # e.g. '90.0', '91.0', '91.1', '_combined'
nc_dir = Ldir['data'] + 'hycom' + tag + '/'

# extrapolation
import time
for vn in vn_list:
    tt0 = time.time()
    print ''
    print 'Filtering ' + vn
    Ofun.process_time_filter(vn, nc_dir)
    print ' Took %.3f seconds' % (time.time() - tt0)
    sys.stdout.flush()   