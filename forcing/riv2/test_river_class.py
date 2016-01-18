"""
Program to test the new river code.
"""

import matplotlib.pyplot as plt
import os
import sys
alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart('cascadia1','base')
import river_class
reload(river_class)

if True:
    # get the list of rivers that we need for a run
    import pandas as pd
    rdf = pd.read_csv(Ldir['run'] + 'rname_list.txt', header=None,
        names=['River Name'])
    rnames = rdf['River Name'].values
    rnames = rnames.tolist()
    # Note, both the Ecology list and the NWS list use green and hamma to
    # refer to these rivers, whereas when we made the rivers for the old
    # MoSSea and PNWTOX models we called them duwamish and hammahamma.
    rnames[rnames.index('duwamish')] = 'green'
    rnames[rnames.index('hammahamma')] = 'hamma'
    # and remove skagit_south
    rnames.pop(rnames.index('skagit_south'))
    # and the fraser is not implemented yet in riv2
    rnames.pop(rnames.index('fraser'))
else:
    # override for testing
    rnames = ['skokomish']

run_type = 'forecast'

if run_type == 'backfill':
    from datetime import datetime
    dt0 = datetime(2015,1,1)
    dt1 = datetime(2015,1,30)
    days = (dt0, dt1)
elif run_type == 'forecast':
    days = ()

plt.close()

for rn in rnames:   
    riv = river_class.River(Ldir)
    riv.name_it(rn)   
    riv.get_ecy_info(riv.name)   
    riv.get_nws_info(riv.name)
    if True:      
        if run_type == 'backfill':
            riv.get_usgs_data(days)
        elif run_type == 'forecast':
            if len(riv.nws_code) > 0:
                riv.get_nws_data()
                if not riv.got_data:
                    riv.get_usgs_data(days)
            else:
                riv.get_usgs_data(days)
    riv.print_info()
   
    if True and not riv.qt.empty:
        fig = plt.figure()
        fig.add_subplot(riv.qt.plot(title=riv.name.title() + ' Flow (' + riv.flow_units + ')', style='-k'))
        plt.show()


