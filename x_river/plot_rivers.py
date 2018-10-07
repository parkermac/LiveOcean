#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot as-run river time series.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc
import argparse
import pandas as pd

import os
import sys
pth = os.path.abspath('../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

Ldir = Lfun.Lstart('cas4', 'v2')

fn = Ldir['LOo'] + 'river/cas4_v2_2017.01.01_2017.12.31.p'

df = pd.read_pickle(fn)