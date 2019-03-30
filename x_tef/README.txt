Code to do TEF extractions and processing.

------------------------------------------------------------------

* ptools/pgrid/tef_section_maker.py is a tool to define sections, if needed.  This is a GUI where you define section endpoints (forced to be either perfectly E-W or N-S), name the section, and define which direction is "landward" (i.e. which way is positive transport).  You can define one or more sections in a session.  At the end, when you push DONE it produces screen output suitable for pasting into new or existing lists in tef_fun.get_sect_df().

------------------------------------------------------------------

* extract_sections.py creates a NetCDF file for each section with arrays of hourly transport and tracer values on the section, arranged as (t, z, x-or-y).  Using command line arguments you can change the run and the day range.  When you run the code you can decide at the command line to extract all the sections from the list, or just one.  The main list is defined in tef_fun.get_sect_df().  Typically this will be done on a remote machine, like boiler, although the defaults are designed to work with model output I have saved on my mac.

Input: ROMS history files over some date range

Output: ptools_output/tef/[*]/extractions/[sect name].nc
where [*] = cas4_v2_lo6biom_2017.07.01_2017.07.03 e.g.
Variables: ocean_time, salt, q, z0, DA0, lon, lat, h, zeta

NOTE: the actual tracer variables to be extracted are defined in tef_fun.start_netcdf.nc().  Not ideal to have it hidden so deep, but it is a convenient place to create the whole list if needed.

------------------------------------------------------------------

* process_sections.py organizes all the transports at each time into salinity bins.

Input: ptools_output/tef/[*]/extractions/[sect name].nc

Output: ptools_output/tef/[*]/processed/[sect name].p
a pickled dict with keys: ['tef_q', 'tef_qs', 'sbins', 'ot', 'qnet', 'fnet', 'ssh']

------------------------------------------------------------------

* bulk_calc.py does the TEF bulk calculations, using the algorithm of Marvin Lorenz, allowing for multiple in- and outflowing layers

Input: ptools_output/tef/[*]/processed/[sect name].p

Output: ptools_output/tef/[*]/bulk/[sect name].p

------------------------------------------------------------------

* bulk_plot.py plots the results of bulk_calc.py, either as a single plot to the screen, or multiple plots to png's.

Input: ptools_output/tef/[*]/bulk/[sect name].p

Output: ptools_output/tef/[*]/bulk_plots/[sect name].p

------------------------------------------------------------------

Code below needs to be updated - not guaranteed to work.

3. Run plot_time_series.py to make plots of time series from each processed station.  You can edit the code to just look at a subset of the stations.

4. Run plot_thalweg_mean.py to make mean values of Qin/Qout, etc. versus along channel distance for a sequence of stations defined by hand in the code.

NOTE: Need to add more lines (e.g. Hood Canal)

Other code:

Run plot_section_map.py to see where all the sections are.  Does not need to be run in order.