#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:13:12 2017

@author: PM5

Code to explore more reliable ways to get the hycom file listings.
"""
exnum = '91.2'
    
from urllib.request import Request, urlopen
from urllib.error import URLError
from socket import timeout
import time
 
xml_name = ('http://tds.hycom.org/thredds/catalog/GLBu0.08/expt_' + 
            exnum + '/forecasts/catalog.xml')  
req = Request(xml_name)
counter = 1
got_file = False
while (counter <= 10) and (got_file == False):
    print('\nCounter = ' + str(counter))    
    tt0 = time.time()
    try:
        xfile = urlopen(req, timeout=30)
    except URLError as e:
        if hasattr(e, 'reason'):
            print(' *We failed to reach a server.')
            print(' -Reason: ', e.reason)
        elif hasattr(e, 'code'):
            print(' *The server couldn\'t fulfill the request.')
            print(' -Error code: ', e.code)
    except timeout:
        print(' *Socket timed out')
    else:
        got_file = True
        print(' Worked fine')
    print(' -took %0.1f seconds' % (time.time() - tt0)) 
    counter += 1

