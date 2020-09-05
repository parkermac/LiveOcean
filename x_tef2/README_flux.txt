README for the flux_ code

This code is focused on the "efflux_reflux" concepts of Cokelet and Stewart (1985, and subsequent).  It picks of the TEF analysis from the point where we have done all the TEF sections, averaged into "bulk" inflows and outflows.  At that point we start working toward a super-simple numerical model in which the tracer concentation in "segments" (the volume in between two or more "sections") can be calculated as dynamically evolving functions of time.

NOTE: (*) means the description/code has not been edited yet for tef2.

------------------------------------------------------------------

* flux_fun.py is a module used throughout this code.

A key piece is the dict flux_fun.segs where we set up the segment names and define (i) the sections on each direction NSEW, and (ii) the rivers flowing into each segment.  This is created by hand while looking at the plots created by flux_seg_map.py below.

A key item is the dict "segs" which whose keys are the segment names, and whose values are in turn dicts which define lists of sections surrounding a segment, and a list of rivers.
segs = {
        'J1':{'S':[], 'N':[], 'W':['jdf1'], 'E':['jdf2'], 'R':['sanjuan', 'hoko']}, etc.
Originally this was created by looking at the plot made by plot_thalweg_mean.py, but the same information is shown more clearly now in the plots made by flux_seg_map.py below.

Also has channel_dict, seg_dict, other things which associate lists of sections or segments with each of four "channels".

flux_fun.make_dist(x,y) makes vectors of distance along lon,lat vectors x,y

------------------------------------------------------------------


------------------------------------------------------------------

* flux_get_vol.py gets the volume (with SSH = 0) of each segment.  Uses a cute "mine sweeper" like algorithm to find all the unmasked rho-grid cells on the grid between the sections that define a segment.

Input: flux_fun.segs, the model grid (hard wired in the code to be cas6), and sect_df = tef_fun.get_sect_df().

Output: LiveOcean_output/tef2/volumes_cas6/...
	volumes.p a pickled DataFrame with segment names for the index and columns: ['volume m3', 'area m2', 'lon', 'lat']
	bathy_dict.p pickled dict of the grid bathy, for plotting
	ji_dict.p pickled dict (keys = segment names) of the ji indices of each segment area, also for plotting
	
(*) flux_V_plot.py makes a simple plot of the segment volumes represented as rectangles:
	LiveOcean_output/tef/misc_figs_cas6/volume_plot.png
	Eventually this might be an interesting framework for plotting the transports used in the flux engine.
	
------------------------------------------------------------------

------------------------------------------------------------------

* flux_get_s.py creates time series of hourly salinity and volume in segments, over a user specified (like a year) period.  Now includes files for the salinity-squared and the net "mixing" (variance destruction due to vertical mixing - need to add horizontal mixing)

Input: ROMS history files

Output: LiveOcean_output/tef2/[*]/flux/hourly_segment_[volume,salinity,mix,salinity2].p, pickled DataFrames where the index is hourly datetime and the columns are the segments.

NOTE: [*] = cas6_v3_lo8b_2017.01.01_2017.12.31, e.g., here and below.

------------------------------------------------------------------

* flux_lowpass_s.py creates time series of daily salinity, volume, and net salt in segments.

Input: LiveOcean_output/tef2/[*]/flux/hourly_segment_[volume,salinity,mix,salinity2].p

Output: LiveOcean_output/tef2/[*]/flux/daily_segment_[volume,salinity,net_salt,mix,salinity2,net_salt2].p

------------------------------------------------------------------

(*) flux_salt_budget.py makes a complete volume and salt budget for a user-specified set of segments.
These budgets are of the form:
	dSnet/dt = QSin + QSout, and
	dV/dt = Qin + Qout + Qr
They are mainly useful for knowing that the budgets add up in each of the basins and years to reasonable accuracy.  The values plotted are DAILY (tidally averaged).

Input: LiveOcean_output/tef2/[*]/flux/daily_segment_[volume,net_salt].p

Output: A screen plot of the budget vs. time, and a saved png such as:
LiveOcean_output/tef2/salt_budget_plots/salt_budget_2017_Salish_Sea.png

------------------------------------------------------------------