# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 06:37:14 2016

@author: PM5

The update backfill code is hangin on this call in hfun.
I want to figure out why.

cd LiveOcean/forcing/hycom

Sample output:

0 sec to get Dataset
127 sec to get 1 3D point
0 sec to get 40 3D points
21 sec to get 4960 3D points
110 sec to get 1240 3D points
0 sec to get 99 3D points
0 sec to get 12276 3D points

=> performance is poor and highly irregular

"""
import netCDF4 as nc
import time
import sys

fn = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1'
#fn = 'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2'

nt0 = 688
nt1 = 718
#nt0 = 10
#nt1 = 40
N = 40
j0 = 1513
j1 = 1637
i0 = 2888
i1 = 2987

tt0 = time.time()
ds = nc.Dataset(fn)
print('%d sec to get Dataset' % (time.time()-tt0))
sys.stdout.flush()
#
#tt0 = time.time()
#x1 = ds.variables['water_temp'][nt0, 0, j0, i0].squeeze()
#print('%d sec to get 1 3D point' % (time.time()-tt0))
#sys.stdout.flush()
#
# tt0 = time.time()
# x2 = ds.variables['water_temp'][nt0, 0:N, j0, i0].squeeze()
# print('%d sec to get %d 3D points' % (time.time()-tt0, x2.size))
# sys.stdout.flush()
#
#tt0 = time.time()
#x3 = ds.variables['water_temp'][nt0, 0:N, j0:j1, i0].squeeze()
#print('%d sec to get %d 3D points' % (time.time()-tt0, x3.size))
#sys.stdout.flush()
#
#tt0 = time.time()
#x4 = ds.variables['water_temp'][nt0:nt1+1, 0:N, j0, i0].squeeze()
#print('%d sec to get %d 3D points' % (time.time()-tt0, x4.size))
#sys.stdout.flush()
#
#tt0 = time.time()
#x5 = ds.variables['water_temp'][nt0, 0, j0, i0:i1].squeeze()
#print('%d sec to get %d 3D points' % (time.time()-tt0, x5.size))
#sys.stdout.flush()
#
# tt0 = time.time()
# x6 = ds.variables['water_temp'][nt0, 0, j0:j1, i0:i1].squeeze()
# print('%d sec to get %d 3D points' % (time.time()-tt0, x6.size))
# sys.stdout.flush()

for nt in range(nt0, nt0+5):#range(nt0, nt1+1):
   tt0 = time.time()
   x7 = ds.variables['water_temp'][nt, 0:N, j0:j1, i0:i1].squeeze()
   print('%d sec to get %d 3D points' % (time.time()-tt0, x7.size))
   sys.stdout.flush()

# and this is the full desired call (for one variable)
# RESULT the call failed
# tt0 = time.time()
# t3d = ds.variables['water_temp'][nt0:nt1+1, 0:N, j0:j1, i0:i1].squeeze()
# print('%d sec to get %d 3D points' % (time.time()-tt0, t3d.size))
# sys.stdout.flush()


ds.close()