"""
Code to print all dt_info files.

"""

# setup
import os
import sys
pth = os.path.abspath('../../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pickle

import hfun
from importlib import reload
reload(hfun)

# initial experiment list
h_list = list(hfun.hy_dict.keys())
h_list.sort()

# specify output directory
out_dir0 = Ldir['data'] + 'hycom2/'
out_dir = out_dir0 + 'dt_lists/'

for hy in h_list:
    print('Experiment = ' + hy)
    f = open(out_dir + 'dt_info_' + hy + '.txt', 'r')
    for line in f:
        sys.stdout.write(' - ' + line)
    f.close()
    print('\n')
