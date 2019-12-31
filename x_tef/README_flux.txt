README for the flux_ code

This code is focused on the "efflux_reflux" concepts of Cokelet and Stewart (1985, and subsequent).  It picks of the TEF analysis from the point where we have done all the TEF sections, averaged into "bulk" inflows and outflows.  At that point we start working toward a super-simple numerical model in which the tracer concentation in "segments" (the volume in between two or more "sections") can be calculated as dynamically evolving functions of time.

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

* flux_make_two_layer.py creates the time-mean two layer versions of the TEF variables, for each section.  These eventually will drive the efflux-reflux integrations.

Input: section "bulk" files like: LiveOcean_output/tef/[*]/bulk/[sect name].nc
where [*] = cas6_v3_lo8b_2017.01.01_2017.12.31 e.g.

Output: LiveOcean_output/tef/[*]/flux/two_layer.p which is a pickled DataFrame whose index is the section names, and whose columns are: ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'].  Here _s and _f indicate that the layer is salty or fresh.  Also we have organized all the fluxes to be positive Eastward or Northward.  The averaging is basically identical to what we did in process_thalweg_mean.py, but here we save the results organized by section instead of the channel lists.

NOTE: Typically the time averaging is done over a whole year.  Eventually we may want to do individual months, but this would require a different organization of the output.

------------------------------------------------------------------

* flux_get_vol.py gets the volume (with SSH = 0) of each segment.  Uses a cute "mine sweeper" like algorithm to find all the unmasked rho-grid cells on the grid between the sections that define a segment.

Input: flux_fun.segs, the model grid (hard wired in the code to be cas6), and sect_df = tef_fun.get_sect_df().

Output: LiveOcean_output/tef/[*]/flux/...
	volumes.p a pickled DataFrame with segment names for the index and columns: ['volume m3', 'area m2', 'lon', 'lat']
	bathy_dict.p pickled dict of the grid bathy, for plotting
	ji_dict.p pickled dict (keys = segment names) of the ji indices of each segment area, also for plotting
	
* flux_V_plot.py makes a simple plot of the segment volumes represented as rectangles:
	LiveOcean_output/tef/[*]/misc_figs/volume_plot.png
	Eventually this might be an interestng framework for plotting the transports used in the flux engine.
	
------------------------------------------------------------------

* flux_seg_map.py

Input: grid, rivers, sections and segments (from volumes.py)

Output: LiveOcean_output/tef/[*]/misc_figs/seg_map_full.png and seg_map_focus.py which are useful for seeing how each segment is related to its defining sections and rivers.  Must be consistent with flux_fun.segs.

------------------------------------------------------------------

* flux_seg_map_simple.py

Input: the output of flux_get_vol.py

Output: LiveOcean_output/tef/[*]/misc_figs/seg_map_simple.png, a graphically appealing plot showing the four basins defined by our "channels" as a colorful map.

------------------------------------------------------------------

* flux_get_A.py is a central and complex part of the flux machinery.  Its job is to take the TEF transports created by flux_make_two_layer.py, and the river flowsm and turn them into efflux-reflux fractions used by the flux engine.

Input: LiveOcean_output/tef/[*]/flux/two_layer.p and...
	river flow extracted over the model year from the forcing (see x_river)
	flux_fun.segs

Output: LiveOcean_output/tef/[*]/flux/q_df.p a pickled DataFrame where the index is the segment list, but with separate salty and fresh entries, like J1_s and J1_f.  The columns are the same, but also there are four more columns at the start: 'ocean_s', 'ocean_f', 'river_s', 'river_f'. From the comments:
# The way q_df is used in the flux_engine is that the row denotes where net transports are
# ending up, and the columns denote where the transports are coming from.  Each entry in
# q_df is a transport (m3/s).  Hence when we multiply a row by the tracer concentrations
# in each segment and then sum along that row, we get the net flux of tracer into or out of
# the segment specified by that row.

NOTE: This assumes that the TEF transports have been averaged over a certain time period, typically a year, and this is done when we make two_layer.p.  But the averaging of river flow is done in this code, so you have to make sure that the averaging periods are the same.

NOTE: Look deep in the code for methods of calculation, and hand adjustments to, the efflux-reflux coefficients.

------------------------------------------------------------------

* flux_w_plot.py makes a plot of vertical velocity in all the segments.

Input: volumes.p and q_df.p created above

Output:LiveOcean_output/tef/[*]/misc_figs/w_plot.png

------------------------------------------------------------------

* *** flux_engine.py ***
This is the main piece of code this whole efflux-reflux analysis has been driving toward.  It takes the transport matrix created by flux_get_A.py, and the volumes, and then uses them to do forward calculations of time-dependent tracer fields.  You can control the boundary conditions (where the tracer comes from), and whether it sinks - like organic particles, using command line arguments.

Input: volumes.p and q_df.p created above

Output: indir = LiveOcean_output/tef/[*]/flux/
	cc.to_pickle(indir + 'cc_' + source + sink_tag + '.p')
	aa.to_pickle(indir + 'aa_' + source + sink_tag + '.p')
	
Here cc is a DataFrame, same index as q_df, and the columns 'c' and 'ca' are the final (after say 6 years) values of the tracer, and an aging tracer.

Also aa is a DataFrame of a full time-dependent array, which you could use to calculate residence times in experiments where you are setting up an initial condition.  For aa the index is time, and the columns are the q_df indices.

------------------------------------------------------------------

* flux_plot_validation.py makes a plot of the steady-state output of the flux engine for the case which is supposed to reproduce the ocean salinity.  This is the only way to test the flux engine and its efflux-reflux coefficients.  It compares the cc_ values to the two_layer.p TEF values.

Input: cc_ocean_salt.p (from flux_engine.py) and two_layer.p

Output: a plot

------------------------------------------------------------------

* flux_plot_cc.py makes a plot of tracer concentration and age from any instance of a cc_ file created by flux_engine.py.

------------------------------------------------------------------

* flux_plot_movie.py makes a movie out of times saved in any instance of an aa_ file created by flux_engine.py.

------------------------------------------------------------------

* flux_get_s.py creates time series of hourly salinity and volume in segments, over a user specified (like a year) period.

Input: ROMS history files

Output: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p, pickled DataFrames where the index is hourly datetime and the columns are the segments.

------------------------------------------------------------------

* flux_lowpass_s.py creates time series of daily salinity, volume, and net salt in segments.

Input: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p

Output: LiveOcean_output/tef/[*]/flux/daily_segment_[salinity,volume,net_salt].p

------------------------------------------------------------------

* flux_salt_budget.py makes a complete salt budget for a user-specified set of segments.

Input: LiveOcean_output/tef/[*]/flux/daily_segment_[net_salt].p

Output: A screen plot of the budget vs. time.

------------------------------------------------------------------
