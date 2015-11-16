"""
Store results of tracker.

Designed for multiple releases.

Result arrays are arrays in (time, particle)

"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path: sys.path.append(alp)
import Lfun; reload(Lfun)
import zfun; reload(zfun) # plotting functions
import cPickle as pickle
import numpy as np

Ldir = Lfun.Lstart()
# NOTE indir should be tracks_2014_CERF for the current experiments
indir = Ldir['LOo'] + 'tracks/'

jdf_list_raw = os.listdir(indir)
jdf_list = []
for m in jdf_list_raw:
    if (m[-2:] == '.p') and 'jdf' in m:
        jdf_list.append(m)
cr_list_raw = os.listdir(indir)
cr_list = []
for m in cr_list_raw:
    if (m[-2:] == '.p') and 'cr' in m:
        cr_list.append(m)
    
# debugging
if False:
    nf = 10
    jdf_list = jdf_list[:nf]
    cr_list = cr_list[:nf]

riv_list = ['jdf', 'cr']
for riv in riv_list:
    if riv == 'jdf':
        fn_list = jdf_list
    elif riv == 'cr':
        fn_list = cr_list    
    # get starting points (assumes these are all reverse runs)
    counter = 0
    NF = len(fn_list)
    pp = dict()
    for fn in fn_list:
        print('%s: working on %d out of %d' % (riv, counter, NF))
        sys.stdout.flush()
        P, G, S, PLdir = pickle.load( open( indir + fn, 'rb' ) )
        if counter == 0:
            # initialize the arrays in the dict
            NP = len(P['lon'][0,:])
            pk = P.keys()
            for vn in pk:
                if vn == 'ot':
                    # first point
                    pp[vn] = P[vn][0]
                    # last point
                    pp[vn + '1'] = P[vn][-1]
                else:
                    pp[vn] = P[vn][0,:].reshape(1,NP)
                    pp[vn + '1'] = P[vn][-1,:].reshape(1,NP)
        else:
            # continue filling the arrays
            for vn in pk:
                if vn == 'ot':
                    pp[vn] = P[vn][0]
                    pp[vn + '1'] = P[vn][-1]
                else:
                    a = pp[vn]
                    b = P[vn][0,:].reshape(1,NP)
                    c = np.concatenate((a,b),axis=0)
                    pp[vn] = c   
                    a1 = pp[vn + '1']
                    b1 = P[vn][-1,:].reshape(1,NP)
                    c1 = np.concatenate((a1,b1),axis=0)
                    pp[vn + '1'] = c1   
        counter += 1
    if riv == 'jdf':
        jdf = pp
    elif riv == 'cr':
        cr = pp    
    
pickle.dump( (jdf, cr, G) , open( indir + 'starters_2014.p', 'wb' ) )

    