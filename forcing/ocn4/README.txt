These programs create the ocean files used to force a roms run, typically for a single day or a 3 day forecast.  They bring in hycom data, either archived extractions for backfill, or from a web source for a forecast.  This is a fairly complicated process, becasue of some lack of regularity of the input formats.

It has been substantially recoded from previous versions to use the new hycom FMRC_best format of the forecast files.

For backfill it will use files of similar format, available from late 2012 to the present, created by the code in hycom2.

===============================================================================
* make_forcing_main.py is the main driver, similar in basic construction and usage to all the code of the same name in other forcing folders (e.g. riv2, tide1, atm).

Input: command line arguments tell it what to do - specifically which gridname_tag = [gtag] to work on, what day [f_string], and whether this is a forecast or backfill.  Then it either gets data from the web, or from LiveOcean_data/hycom1.

Output: LiveOcean_output/[gtag]/[f_string]/ocn4/...
(1) Info/process_status.csv and screen_out.txt (great for debugging)
(2) Data/...
- data.nc (all the hycom data for this period - only for the new FMRC_best forecast)
- coord_dict.p (pickled dict of coordinate vectors - just like in hycom1)
- h[data_string].p (pickled dict of fields for a day - for a forecast there will be 8 of these)
- fh[date_string].p (time filterd version of h[].p files - for a forecast there will be 4 of these)
- xfh[date_string].p (spatially extrapolated version of fh[].p files - for a forecast there will be 4 of these, no nans)
(3) ocean_[clm, ini, bry].nc are the files for forcing roms

Notes:
- If something goes wrong with getting the data for a forecast, then it uses "planB" in which the ocean_clm.nc file is copied from the day before, and one day is added to its last time.
- in backfill mode it uses the archives files , e.g., LiveOcean_data/hycom2/hy6/h2019.01.01.nc.  Currently no planB for this operation.

Modules:
- Ofun.py: the main workhorse functions to get data, filter in time, and extrapolate
- Ofun_bio.py: add and fill bio variables, e.g. using regressions against salinity
- Ofun_CTD.py: fill fields in the Salish Sea and coastal estuaries using CTD observations.  Currently only set to work for 1/1/2017.  The goal is to get a better initial condition.
- Ofun_nc.py: functions for making the ouput NetCDF files (3).

===============================================================================
* check_ocn.py plot ocean fields to compare raw hycom, extrapolated hycom, and interpolated roms.

Input: a selection of the files for a given [f_string]/ocn3, as well as the grid.nc file from LiveOcean_data/grids.

Output: a screen plot of map fields for a variable - useful for debugging

===============================================================================
* check_bio.py plot ocean fields to check on interpolated roms fields of the bio variables.

Input: ocean_ini.nc a given [f_string]/ocn3, as well as the grid.nc file from LiveOcean_data/grids.

Output: a screen plot of map fields for a bunch of variables - useful for debugging