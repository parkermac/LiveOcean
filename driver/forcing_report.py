"""
This generates a report on the results of all the forcing-generation code.
"""

# setup
import os; import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

f_list = ['atm', 'ocn', 'riv', 'tide']
#f_list = ['tide']

for which_force in f_list:
    Info = Lfun.csv_to_dict(Ldir['LOo'] + 'current_Info/'
    + which_force + '/Info_for_worker.csv')
    #ifw = Lfun.csv_to_dict(Info['Fdir_fi'] + 'Info_from_worker.csv')
    ifw = Lfun.csv_to_dict(Ldir['LOo'] + Ldir['gtag'] + '/'
    + Info['f_string'] + '/' + which_force + '/Info/' 'Info_from_worker.csv')
    
    print '/n******************************************************'
    template = '%20s:%30s'
    print template % ('results from', which_force)
    i_list = ['result', 'start_time', 'end_time', 'var_start_time', 'var_end_time']
    for ii in i_list:
        try:
            print template % (ii, ifw[ii])
        except:
            print template % (ii, 'MISSING')


