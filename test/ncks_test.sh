#!/bin/bash

#time ncks -D8 -d time,2020-06-21T00:00,2020-06-29T00:00,8 -d lon,229.,239.,1 -d lat,39.,53.,1 -v surf_el,water_temp,salinity,water_u,water_v https://tds.hycom.org/thredds/dodsC/GLBy0.08/latest -4 -O /Users/pm8/Desktop/ncks_test.nc


time ncks -d time,2020-06-21T00:00,2020-06-28T00:00,8 -d lon,229.,239.,1 -d lat,39.,53.,1 -v surf_el,water_temp,salinity,water_u,water_v https://tds.hycom.org/thredds/dodsC/GLBy0.08/latest -4 -O /Users/pm8/Desktop/ncks_test.nc


# GLBy0.08_latest_20200609_20200615_`date +%s`.nc4

# time ncks -D8 -d time,2020-06-22T00:00,2020-06-25T00:00,8 -d lon,229.,239.,1 -d lat,39.,53.,1 -v surf_el,water_temp,salinity,water_u,water_v https://beta.hycom.org/thredds/dodsC/GLBy0.08/latest -4 -O GLBy0.08_latest_20200609_20200615_`date +%s`.nc4