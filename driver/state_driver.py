#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 07:41:32 2017

@author: PM5

This is a driver for the LiveOcean daily forecast system.  It is meant
to be run many times per day by cron, on fjord or gaggle, in order to
run the forecast in the most efficient and flexible way possible.

"""

import os
import sys
from datetime import datetime, timedelta
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun

import subprocess

time_format = '%Y.%m.%d %H:%M:%S'

# set defaults
gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio1'
run_type = 'forecast'  # backfill or forecast
start_type = 'continuation'  # new or continuation
date_string = datetime.now().strftime(format='%Y.%m.%d')

Ldir = Lfun.Lstart(gridname, tag)
Ldir['date_string'] = date_string
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name
Ldir['run_type'] = run_type
Ldir['start_type'] = start_type
Ldir['date_string'] = date_string
Ldir['ex_name'] = ex_name

# Make the directory tree for forcing.
Ldir['LOog'] = (Ldir['LOo'] + Ldir['gtag'] + '/')
Ldir['LOogf'] = (Ldir['LOog'] + 'f' + date_string + '/')
# Make the directory tree for model output.
Ldir['LOro'] = (Ldir['roms'] + 'output/')
Ldir['LOrog'] = (Ldir['LOro'] + Ldir['gtagex'] + '/')
Ldir['LOrogf'] = (Ldir['LOrog'] + 'f' + date_string + '/')

# Make directory for state files
states_dir = Ldir['LO'] + 'driver/states/'
Lfun.make_dir(states_dir)

# create a dict of tasks and their requirements
task_dict = {'tide':[],
            'riv':[],
            'ocn':[],
            'bio':['ocn'],
            'atm':[],
            'roms':['riv', 'tide', 'bio', 'atm'],
            'carbon':['roms'],
            'low_pass':['carbon'],
            'azu':['carbon']}

# create the full path to the state.csv file for this day
state_fn = (states_dir + 'state_' + Ldir['date_string'] + '.csv')

# initialize state.csv if it does not exist
try:
    state_dict = Lfun.csv_to_dict(state_fn)
except FileNotFoundError:
    state_dict = dict()
    for task in task_dict.keys():
        state_dict[task] = 'to_do'    
    Lfun.dict_to_csv(state_dict, state_fn)

# step through the tasks
for task in task_dict.keys():
    
    # get the state
    state_dict = Lfun.csv_to_dict(state_fn)
    
    if state_dict[task] == 'to_do':
        req_all_done = [True]
        requirement_list = task_dict[task]
        for requirement in requirement_list:
            if state_dict[requirement] == 'done':
                pass
            elif state_dict[requirement] == 'to_do':
                req_all_done.append(False)
        if req_all_done.count(False) == 0:
            # run this task
            print('Running ' + task)
            if task == 'tide':
                command_list = [Ldir['LO']+'driver/driver_forcing1.sh',
                            '-g', Ldir['gridname'], '-t', Ldir['tag'],
                            '-f', task, '-r', Ldir['run_type'], '&']
                proc = subprocess.Popen(command_list,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = proc.communicate()
                if False:
                    print('Screen output:')
                    print(out.decode())
                    print('Error messages:')
                    print(err.decode())
                
                # check on the result and update the state if it worked
                status_dict = Lfun.csv_to_dict(Ldir['LOogf'] + task
                            + '/Info/process_status.csv')
                if status_dict['result'] == 'success':
                    print('updating state of task ' + task)
                    state_dict[task] = 'done'
                    Lfun.dict_to_csv(state_dict, state_fn)
                    
    elif state_dict[task] == 'done':
        pass
