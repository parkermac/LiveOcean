"""
Code to plot a dt_log_ file.
"""

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

import pandas as pd
import matplotlib.pyplot as plt

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'extract/'

# choose the mooring extraction to plot
print('\n%s\n' % '** Choose dt_log file to plot **')
m_list_raw = os.listdir(indir)
m_list_raw.sort()
m_list = [m for m in m_list_raw if (('.txt' in m) and ('dt_log_' in m))]
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
log_file = m_dict[my_npt]
fn = indir + log_file

df = pd.read_csv(fn,parse_dates=True,index_col='Date')

plt.close('all')

df.plot(title=log_file.strip('.txt'))

plt.show()
