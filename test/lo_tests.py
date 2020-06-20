"""
Module of code to test various aspects of the LiveOcean system.
"""

# imports
import os, sys
sys.path.append(os.path.abspath('../alpha'))
import Lfun

# things we needtest
Ldir = Lfun.Lstart()

def check_Ldir():
    print('Check Ldir'.center(50,'='))
    for k in Ldir.keys():
        print('%15s : %s' % (k, Ldir[k]))
        
def check_read_his():
    print('Check read his'.center(50,'='))
    fn = Ldir['data'] + 'test/ocean_his_0001.nc'
    try:
        import netCDF4 as nc
        ds = nc.Dataset(fn)
        print(ds['salt'])
    except Exception as e:
        print('Exception:')
        print(e)
        
def check_azure():
    print('Check azure'.center(50,'='))
    try:
        os.system('python ' + Ldir['LO'] + 'x_misc/move_to_azure.py')
    except Exception as e:
        print('Exception:')
        print(e)
    
    


if __name__ == '__main__':
    # this will run when you just run the module
    
    check_Ldir()
    
    check_read_his()
    
    check_azure()
