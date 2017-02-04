import Ofun
from importlib import reload
reload(Ofun)

# create a list of url's of the preprocessed HYCOM files for this forecast
cc = 0
cc_max = 10
while cc < cc_max:
    try:
        print(' attempt number ' + str(cc))
        fn_list = Ofun.get_hycom_file_list()
        print('** DONE getting catalog')
        cc = cc_max
    except:
        cc += 1
# get a selection of the raw list (e.g. one per day)
varf_dict, dt_list2 = Ofun.get_varf_dict(fn_list)
var_list = list(varf_dict.keys())
NT = len(dt_list2)



    
    
        