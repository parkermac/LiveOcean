This code is designed to make a backfill archive for the GOFS 3.1 version of HYCOM:

The newer GOFS 3.1  goes back as Renanalysis to 1994, and as Analysis to the present.  What we use in forcing/ocn3 is "GLBy0.08 grid is 0.08 deg lon x 0.04 deg lat that covers 80S to 90N." so this is a different grid than GLBu, having twice the latitude resolution.

Start here to see how to format request strings:
https://ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0/dataset.html

=====================================================================


=====================================================================
**** ALL FILES BELOW ARE OLD NOTES, FROM hycom1 *********************

These files are for pulling in and processing backfill HYCOM data. Note that they are currently configured only to pull in the GLBu version of the output, which was turned off after 2018.11.20.  I have been running this on fjord so far, because it is the home for LiveOcean_data.  I guess in principle it would work from boiler as well, because it looks in the same place as fjord for _data.

Since GLBu is no longer being added to, there should be no need to re-run the code in this folder (2018.12.25).

NOTE: need to figure out how to deal with backfill for the newer GLBy format.

=====================================================================
* get_hycom_days.py it the main driver for filling out the archive.  While it can get all the daily extractions for all the experiment numbers, in practice it is much more parsimonious, only getting an extraction if it is not already in the archive.

Input: it goes out to a url like:
http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_ + ['90.9', '91.0', '91.1', '91.2']
where the things like 91.2 are the various hycom experiments, each of which fills out a certain span of time with one extraction per day.

Output:
LiveOcean_data/hycom1/h2016.12.25.p (for example; one for each day) a pickled dict containing:
['v3d', 't3d', 'result', 'ssh', 'u3d', 's3d', 'dt'] where result is something like 'success' and dt is a datetime (at zero hour of the day, presumably UTC).  The dynamical variables are 2D or 3D arrays.  Importantly, they are MASKED ARRAYS and a lot of how I treat them later in ocn3 assumed this to be the case.  A problem that arose was that in the newer forecast fields (GLBy that I get using ncss - see forcing/ocn3) they show up as regular nd arrays with nans, not masked arrays.  This caused a lot of trouble.  All are packed in my standard way [z, y, x] or [y, x].  Velocities are m/s, ssh is m, salinity is psu, and temperature, t3d, is deg C in-situ temperature (we convert to potential temperature later, e.g. in forcing/ocn3, because this is what ROMS uses).

It takes about 30 sec per day on fjord.

=====================================================================
* hycom_listing.py gives a summary of what is in the hycom archive for all exnums.  It also creates the coords_dict.p.

Input: it goes out to a url like:
http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_ + ['90.9', '91.0', '91.1', '91.2']

Output:
(1)LiveOcean_data/hycom1/dt_dict.p which is a dict of lists of the datetimes of the times (one per day) that exist for each exnum.
(2) LiveOcean_data/hycom1/coords_dict.p which is a pickled dict of dicts, e.g. coord_dict = coords_dict['91.2']
	and coord_dict contains:
	['j1', 'z', 'i1', 'j0', 'lat', 'i0', 'lon'] which are either 1D arrays lon, lat, z , or integers (i0, etc.) that give the location we extracted from the hycom grid.  As it happens all four coord_dicts are the same, so you can always safely use the last one.  Note that z is 40 levels from -5000 to 0 with spacing anything from 1000 m (deep est 2 layers) to 2 m (top 6 layers). lat has 150 values, and lon has 100.
(3)And screen output containing:
Working on 90.9
 dt start = 2012-01-25 00:00:00
 dt end   = 2013-08-20 00:00:00
 Have 440 out of 573 days (missing 133)
 Indices 2887:2987 1487:1637

Working on 91.0
 dt start = 2013-08-17 00:00:00
 dt end   = 2014-04-08 00:00:00
 Have 234 out of 234 days (missing 0)
 Indices 2887:2987 1487:1637

Working on 91.1
 dt start = 2014-04-07 00:00:00
 dt end   = 2016-04-18 00:00:00
 Have 738 out of 742 days (missing 4)
 Indices 2887:2987 1487:1637

Working on 91.2
 dt start = 2016-04-18 00:00:00
 dt end   = 2018-11-20 00:00:00
 Have 920 out of 946 days (missing 26)
 Indices 2887:2987 1487:1637
 
It runs in just a few seconds on fjord.
 
===================================================================== 
NOTES on HYCOM versions:
 
GLBu0.08 is a version of the hycom output that has been interpolated to a regular lon, lat grid and regular vertical levels.
From the website: https://www.hycom.org/data/glbu0pt08
"Native hycom .[ab] data converted to NetCDF interpolated to a uniform 0.08 degree lat/lon grid between 80.48S and 80.48N and interpolated to 40 standard z-levels. Five fields are provided: SSH, eastward velocity, northward velocity, in-situ temperature, and salinity."

Looking here:
https://www.hycom.org/dataserver you can see the time ranges of the various products.

This code was developed for the GOFS 3.0 output, which goes back as Reanalysis to 1992, and extends as Analysis to 2018.11.20.

Then the newer GOFS 3.1  goes back as Renanalysis to 1994, and as Analysis to the present.  What we use in forcing/ocn3 is "GLBy0.08 grid is 0.08 deg lon x 0.04 deg lat that covers 80S to 90N." so this is a different grid than GLBu, having twice the latitude resolution.
