"""
Horizontal extrapolation of HYCOM NetCDF files.

Do this BEFORE filtering.
"""
# setup
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

ofp = os.path.abspath('../ocn')
if ofp not in sys.path:
    sys.path.append(ofp)
import Ofun

vn_list = ['ssh', 't3d', 's3d', 'u3d', 'v3d']
tag = '91.1'
nc_dir = Ldir['data'] + 'hycom' + tag + '/'

# extrapolation
import time
for vn in vn_list:
    tt0 = time.time()
    print('\nExtrapolating ' + vn)
    sys.stdout.flush()
    Ofun.process_extrap(vn, nc_dir)
    print(' Took %.3f seconds' % (time.time() - tt0))
    sys.stdout.flush()    