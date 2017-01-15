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
var_list = ['ssh','ts3z','uv3z']    
varf_dict = dict()
for var_name in var_list:
    varf_dict[var_name] = Ofun.make_shortened_list(fn_list, var_name)
# check that all the lists are the same length
list_len = []
for var_name in var_list:
    list_len.append(len(varf_dict[var_name]))
if list_len.count(list_len[0]) != len(list_len):
    print('WARNING: Lists in varf_dict are different lengths!')
    NT = min(list_len) # an attempt to dea with the situation when
    # some hycome variables have fewer data records than others
else:
    NT = list_len[0] # the number of times, to use below   
