Code to do TEF extractions and processing.

------------------------------------------------------------------

* ptools/pgrid/tef_section_maker.py is a tool to define sections, if needed.  This is a GUI where you define section endpoints (forced to be either perfectly E-W or N-S), name the section, and define which direction is "landward" (i.e. which way is positive transport).  You can define one or more sections in a session.  At the end, when you push DONE it produces screen output suitable for pasting into new or existing lists in tef_fun.get_sect_df().

------------------------------------------------------------------
NOTE: all code below is in the directory LiveOcean/x_tef/
------------------------------------------------------------------

* tef_fun.py module for this code.  Includes tef_fun.get_sect_df() which returns a DataFrame of all the section names and their lon,lat endpoints.  The entries in this are pasted in by hand while using tef_section_maker.py as documented above.  Then to control which sections will be processed by later programs (like to do an extraction) you comment/uncomment lines in this function.  This is not great coding - it might be better to have one function that returns the info for any section name, and then make specific lists for specific projects.

------------------------------------------------------------------

* extract_sections.py creates a NetCDF file for each section with arrays of hourly transport and tracer values on the section, arranged as (t, z, x-or-y).  Using command line arguments you can change the run and the day range.  When you run the code you can decide at the command line to extract all the sections from the list, or just one.  The main list is defined in tef_fun.get_sect_df().  Typically this will be done on a remote machine, like boiler, although the defaults are designed to work with model output I have saved on my mac.

Input: ROMS history files over some date range

Output: LiveOcean_output/tef/[*]/extractions/[sect name].nc
where [*] = cas4_v2_lo6biom_2017.07.01_2017.07.03 e.g.
Variables: ocean_time, salt, q, z0, DA0, lon, lat, h, zeta

NOTE: the actual tracer variables to be extracted are defined in tef_fun.start_netcdf.nc().  Not ideal to have it hidden so deep, but it is a convenient place to create the whole list if needed.

NOTE: the sign conventions are defined by the sign of the column "landward" in sect_df.  This is a hangover from when I was trying to be clever about the direction of inflow.  In restospect it would have been better to just use eastward or northward positive.

------------------------------------------------------------------

* process_sections.py organizes all the transports at each time into salinity bins.

Input: LiveOcean_output/tef/[*]/extractions/[sect name].nc

Output: LiveOcean_output/tef/[*]/processed/[sect name].p
a pickled dict with keys: ['tef_q', 'tef_qs', 'sbins', 'ot', 'qnet', 'fnet', 'ssh']
I think these are defined as:
	tef_q transport in salinity bins, hourly, (m3/s)
	tef_qs salt transport in salinity bins, hourly, (m3/s)
	sbins the salinity bin centers
	ot ocean time (sec from 1/1/1970)
	qnet section integrated transport (m3/s)
	fnet section integrated tidal energy flux (units?)
	ssh section averaged ssh (m)
	
------------------------------------------------------------------

* bulk_calc.py does the TEF bulk calculations, using the algorithm of Marvin Lorenz, allowing for multiple in- and outflowing layers

Input: LiveOcean_output/tef/[*]/processed/[sect name].p

Output: LiveOcean_output/tef/[*]/bulk/[sect name].p

------------------------------------------------------------------

* bulk_plot.py plots the results of bulk_calc.py, either as a single plot to the screen, or multiple plots to png's.

Input: LiveOcean_output/tef/[*]/bulk/[sect name].p

Output: LiveOcean_output/tef/[*]/bulk_plots/[sect name].p

------------------------------------------------------------------

* plot_physical_section.py plots a user-selected section as time-averaged properties in x/y-z space, just like a standard Eulerian average.

------------------------------------------------------------------

------------------------------------------------------------------

------------------------------------------------------------------



Code below needs to be updated - not guaranteed to work.

3. Run plot_time_series.py to make plots of time series from each processed station.  You can edit the code to just look at a subset of the stations.

4. Run plot_thalweg_mean.py to make mean values of Qin/Qout, etc. versus along channel distance for a sequence of stations defined by hand in the code.

NOTE: Need to add more lines (e.g. Hood Canal)

Other code:

Run plot_section_map.py to see where all the sections are.  Does not need to be run in order.