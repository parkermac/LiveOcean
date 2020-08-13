"""
Module of code to test various aspects of the LiveOcean system.
"""

# imports
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun

Ldir = Lfun.Lstart()

def print_head(check_str):
    print('')
    print(('Check '+check_str).center(70,'='))
    print('')

def check_Ldir():
    print_head('Ldir')
    for k in Ldir.keys():
        print('%15s : %s' % (k, Ldir[k]))
        
def check_read_his():
    print_head('NetCDF')
    fn = Ldir['data'] + 'test/ocean_his_0001.nc'
    try:
        import netCDF4 as nc
        ds = nc.Dataset(fn)
        print(ds['salt'])
    except Exception as e:
        print('Exception:')
        print(e)
        
def check_azure():
    print_head('Azure')
    try:
        os.system('python ' + Ldir['LO'] + 'x_misc/move_to_azure.py')
    except Exception as e:
        print('Exception:')
        print(e)
        
def check_utide():
    
    print_head('UTide')
    try:
        import utide
        import pandas as pd
        from matplotlib.dates import date2num
        df = pd.read_pickle(Ldir['data'] + 'test/tide_9447130_2017.p')
        t = date2num(df.index.to_pydatetime())
        z = df['eta'].to_numpy()
        h = utide.solve(t, z, v=None,
                     lat=47.6,
                     nodal=False,
                     trend=False,
                     method='ols',
                     conf_int='linear',
                     Rayleigh_min=0.95)
        f = h.aux.frq[h.name == 'M2'][0]
        print('M2 period = ' + str(1/f) + ' hours')
        # h.aux.freq has units cyles/hour
        # so for f = h.aux.frq[h.name == 'M2'][0] we get
        # 1/f = 12.420601202671868 (hours per cycle)
        # h.A is amplitude (m), h.g is phase (degrees)
        
    except Exception as e:
        print('Exception:')
        print(e)

    
if __name__ == '__main__':
    # this will run when you just run the module
    
    check_Ldir()

    check_read_his()

    check_azure()
    
    check_utide()
