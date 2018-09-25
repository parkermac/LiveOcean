README for the TEF code, Parker MacCready, 2018.09.19 and after

Steps and notes

1. Run extract_sections.py, which creates a NetCDF file for each section with arrays of hourly transport and tracer values on the section, arranged as (t, z, x-or-y).  Using command line arguments you can change the run and the day range.  You can edit the code by hand to just extract a subset of the sections.  The main list is defined in tef_fun.get_sect_df().  Typically this will be done on a remote machine, like boiler, although the defaults are designed to work with model output I have saved on my mac.

NOTE: the actual tracer variables to be extracted are defined in tef_fun.start_netcdf.nc().  Not ideal to have it hidden so deep, but it is a convenient place to create the whole list if needed.


2. Run process_sections.py, which organizes all the transports at each time into salinity bins.

3. Run plot_time_series.py to make plots of time series from each processed station.  You can edit the code to just look at a subset of the stations.

4. Run plot_thalweg_mean.py to make mean values of Qin/Qout, etc. versus along channel distance for a sequence of stations defined by hand in the code.

NOTE: Need to add more lines (e.g. Hood Canal)

Other code:

Run plot_section_map.py to see where all the sections are.  Does not need to be run in order.