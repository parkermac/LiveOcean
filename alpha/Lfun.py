"""
Functions for LiveOcean.
"""

def Lstart(gridname='BLANK', tag='BLANK'):
    """
    This is to set environment variables in the LiveOcean system
    using values in a csv file.  It is similar to Lstart.m in the
    MATLAB code, but it returns a dictionary instead of a structure.

    We use input parameters to allow for different gridnames and tags to
    coexist.
    """
    # put top level information from input into a dict
    Ldir = dict()
    Ldir['gridname'] = gridname
    Ldir['tag'] = tag

    # Build information on the directory structure.
    import os
    which_home = os.environ.get("HOME") # This works even when called by cron.
    which_host = os.environ.get('HOSTNAME')
    
    if 'Parkers-MacBook-Pro' in which_host: # mac version
        Ldir['env'] = 'pm_mac'
        Ldir['parent'] = which_home + '/Documents/'
        Ldir['roms'] = Ldir['parent'] + 'LiveOcean_roms/'
        Ldir['which_matlab'] = '/Applications/MATLAB_R2017a.app/bin/matlab'
        
    elif (which_home == '/home/parker') and ('fjord' in which_host):
        Ldir['env'] = 'pm_fjord'
        Ldir['parent'] = '/data1/parker/'
        Ldir['roms'] = '/pmr1/parker/LiveOcean_roms/'
        Ldir['which_matlab'] = '/usr/local/bin/matlab'
        
    elif (which_home == '/home/parker') and ('boiler' in which_host):
        Ldir['env'] = 'pm_boiler'
        Ldir['parent'] = '/data1/parker/'
        Ldir['roms'] = '/pmr1/parker/LiveOcean_roms/'
        Ldir['which_matlab'] = '/usr/local/bin/matlab'
        
    elif which_home == '/usr/lusers/darrd':
        Ldir['env'] = 'dd_mox'
        Ldir['parent'] = '/gscratch/macc/darrd/LOcean2/'
        Ldir['roms'] = '/gscratch/macc/darrd/LOcean2/LiveOcean_roms/'
        
    elif which_home == '/usr/lusers/pmacc':
        Ldir['env'] = 'pm_mox'
        Ldir['parent'] = '/gscratch/macc/parker/'
        Ldir['roms'] = '/gscratch/macc/parker/LiveOcean_roms/'
        
    else:
        print('Trouble filling out environment variables in Ldir')

    # and add a few more things
    Ldir['gtag'] = Ldir['gridname'] + '_' + Ldir['tag']
    Ldir['LO'] = Ldir['parent'] + 'LiveOcean/'
    Ldir['LOo'] = Ldir['parent'] + 'LiveOcean_output/'
    Ldir['data'] = Ldir['parent'] + 'LiveOcean_data/'
    Ldir['grid'] = Ldir['data'] + 'grids/' + Ldir['gridname'] + '/'
    Ldir['forecast_days'] = 3

    return Ldir

def make_dir(dirname, clean=False):
    # Make a directory if it does not exist.
    # Use clean=True to clobber the existing directory.
    import os
    if clean == True:
        import shutil
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass # assume OSError was raised because directory already exists

def dict_to_csv(dict_name, csv_name, write_mode='w'):
    # Write the contents of a dict to a two-column csv file.
    # The write_mode can be wb (overwrite, binary) or ab (append binary).
    # Binary mode is better across platforms.
    # 2015.11.25 changed to 'w' because it was throwing an error in python 3
    import csv
    with open(csv_name, write_mode) as ff:
        ww = csv.writer(ff)
        for kk in dict_name.keys():
            ww.writerow((kk, dict_name[kk]))

def csv_to_dict(csv_name):
    # Reads two-column csv file into a dict.
    import csv
    dict_name = dict()
    with open(csv_name) as ff:
        for row in csv.reader(ff):
            dict_name[row[0]] = row[1]
    return dict_name

def run_worker(Ldir, worker_type='matlab'):
    # run the worker code using subprocess
    if worker_type == 'matlab':
        # pass arguments to a matlab program
        import subprocess
        func = ("make_forcing_worker(\'" +
            Ldir['gridname'] + "\',\'" +
            Ldir['tag'] + "\',\'" +
            Ldir['date_string'] + "\',\'" +
            Ldir['run_type'] + "\',\'" +
            Ldir['LOogf_f'] + "\')")
        cmd = Ldir['which_matlab']
        run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]        
        proc = subprocess.run(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('\n-main: screen output from worker-')
        print(proc.stdout.decode())
    else:
        print('other worker types not implemented yet')


def datetime_to_modtime(dt):
    # This is where we define how time will be treated
    # in all the model forcing files.
    # NOTE: still need to define time zone (UTC?)

    # INPUT:
    # dt is a single datetime value

    # this returns seconds since 1/1/1970
    from datetime import datetime
    t = (dt - datetime(1970,1,1,0,0)).total_seconds()
    return t

def modtime_to_datetime(t):
    # input seconds since 1/1/1970 (single number)
    from datetime import datetime, timedelta
    dt = datetime(1970,1,1,0,0) + timedelta(seconds=t)
    return dt

def modtime_to_mdate_vec(mt_vec):
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates

    # input numpy vector of seconds since 1/1/1970

    # first make a list of datetimes
    dt_list = []
    for mt in mt_vec:
        dt_list.append(datetime(1970,1,1,0,0) + timedelta(seconds=mt))

    md_vec = mdates.date2num(dt_list)

    return md_vec


