"""
Shared helper functions for the forcing code.
"""
def intro():
    # setup
    import os
    import sys
    import argparse
    from datetime import datetime
    # This relative path to alpha is meant to work only when intro()
    # is called from the forcing directories (not great coding).
    alp = os.path.abspath('../../alpha')
    if alp not in sys.path:
        sys.path.append(alp)
    # Note: the path "alp" will now also work for the calling function
    import Lfun

    # set defaults
    gridname = 'cascadia1'
    tag = 'base'
    cwd = os.getcwd()
    icwd = cwd.rfind('/')
    frc = cwd[icwd+1:]
    run_type = 'forecast'  # backfill or forecast
    start_type = 'continuation'  # new or continuation
    # Example of date_string is 2015.09.19
    date_string = datetime.now().strftime(format='%Y.%m.%d')
    ex_name = 'lo1'

    # optional command line arguments, can be input in any order
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gridname', nargs='?', const=gridname, type=str, default=gridname)
    parser.add_argument('-t', '--tag', nargs='?', const=tag, type=str, default=tag)
    parser.add_argument('-f', '--frc', nargs='?', const=frc, type=str, default=frc)
    parser.add_argument('-r', '--run_type', nargs='?', const=run_type, type=str, default=run_type)
    parser.add_argument('-s', '--start_type', nargs='?', const=start_type, type=str, default=start_type)
    parser.add_argument('-d', '--date_string', nargs='?', const=date_string, type=str, default=date_string)
    parser.add_argument('-x', '--ex_name', nargs='?', const=ex_name, type=str, default=ex_name)
    args = parser.parse_args()

    # get the dict Ldir
    Ldir = Lfun.Lstart(args.gridname, args.tag)
    Ldir['date_string'] = args.date_string
    Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

    # add the arguments to Ldir, because some code needs them
    Ldir['gridname'] = args.gridname
    Ldir['tag'] = args.tag
    Ldir['frc'] = args.frc
    Ldir['run_type'] = args.run_type
    Ldir['start_type'] = args.start_type
    Ldir['date_string'] = args.date_string
    Ldir['ex_name'] = args.ex_name

    # Make the directory tree for this forcing, if needed. This is redundant
    # with what the driver does (except that it clobbers nothing), and is
    # only included here so that we can test the python code without using
    # the driver.
    Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
    Ldir['LOogf'] = (Ldir['LOog'] + 'f' + args.date_string + '/')
    Ldir['LOogf_f'] = (Ldir['LOogf'] + args.frc + '/')
    Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
    Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')
    Lfun.make_dir(Ldir['LOog'])
    Lfun.make_dir(Ldir['LOogf'])
    Lfun.make_dir(Ldir['LOogf_f'])
    Lfun.make_dir(Ldir['LOogf_fi'])
    Lfun.make_dir(Ldir['LOogf_fd'])

    # screen output
    print('MAIN: frc = ' + args.frc + ', run_type = ' + args.run_type
        + ', date_string = ' + args.date_string)
    print('MAIN start time = ' + str(datetime.now()))

    return Ldir, Lfun

def intro_2():
    # The purpose of this version is just to handle the new placement of the
    # dot_in forcing files (inside directories named for grids). 8/2/2016
    #
    # setup
    import os
    import sys
    import argparse
    from datetime import datetime
    # This relative path to alpha is meant to work only when intro()
    # is called from one below the forcing directories (REALLY not great coding).
    alp = os.path.abspath('../../../alpha')
    if alp not in sys.path:
        sys.path.append(alp)
    # Note: the path "alp" will now also work for the calling function
    import Lfun

    # set defaults
    gridname = 'cascadia1'
    tag = 'base'
    cwd = os.getcwd()
    icwd = cwd.rfind('/')
    frc = cwd[icwd+1:]
    run_type = 'forecast'  # backfill or forecast
    start_type = 'continuation'  # new or continuation
    # Example of date_string is 2015.09.19
    date_string = datetime.now().strftime(format='%Y.%m.%d')
    ex_name = 'lo1'

    # optional command line arguments, can be input in any order
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gridname', nargs='?', const=gridname, type=str, default=gridname)
    parser.add_argument('-t', '--tag', nargs='?', const=tag, type=str, default=tag)
    parser.add_argument('-f', '--frc', nargs='?', const=frc, type=str, default=frc)
    parser.add_argument('-r', '--run_type', nargs='?', const=run_type, type=str, default=run_type)
    parser.add_argument('-s', '--start_type', nargs='?', const=start_type, type=str, default=start_type)
    parser.add_argument('-d', '--date_string', nargs='?', const=date_string, type=str, default=date_string)
    parser.add_argument('-x', '--ex_name', nargs='?', const=ex_name, type=str, default=ex_name)
    args = parser.parse_args()

    # get the dict Ldir
    Ldir = Lfun.Lstart(args.gridname, args.tag)
    Ldir['date_string'] = args.date_string
    Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

    # add the arguments to Ldir, because some code needs them
    Ldir['gridname'] = args.gridname
    Ldir['tag'] = args.tag
    Ldir['frc'] = args.frc
    Ldir['run_type'] = args.run_type
    Ldir['start_type'] = args.start_type
    Ldir['date_string'] = args.date_string
    Ldir['ex_name'] = args.ex_name

    # Make the directory tree for this forcing, if needed. This is redundant
    # with what the driver does (except that it clobbers nothing), and is
    # only included here so that we can test the python code without using
    # the driver.
    Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
    Ldir['LOogf'] = (Ldir['LOog'] + 'f' + args.date_string + '/')
    Ldir['LOogf_f'] = (Ldir['LOogf'] + args.frc + '/')
    Ldir['LOogf_fi'] = (Ldir['LOogf_f'] + 'Info/')
    Ldir['LOogf_fd'] = (Ldir['LOogf_f'] + 'Data/')
    Lfun.make_dir(Ldir['LOog'])
    Lfun.make_dir(Ldir['LOogf'])
    Lfun.make_dir(Ldir['LOogf_f'])
    Lfun.make_dir(Ldir['LOogf_fi'])
    Lfun.make_dir(Ldir['LOogf_fd'])

    # screen output
    print('MAIN: frc = ' + args.frc + ', run_type = ' + args.run_type
        + ', date_string = ' + args.date_string)
    print('MAIN start time = ' + str(datetime.now()))

    return Ldir, Lfun

def finale(result_dict, Ldir, Lfun):

    # write results to an output file for the driver
    csv_name_out = Ldir['LOogf_fi'] + 'process_status.csv'
    Lfun.dict_to_csv(result_dict, csv_name_out)

    from datetime import datetime
    print('MAIN end time = ' + str(datetime.now()))

