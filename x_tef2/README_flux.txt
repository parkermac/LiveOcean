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

(*) flux_make_two_layer.py creates the time-mean two layer versions of the TEF variables, for each section.  These eventually will drive the efflux-reflux integrations.

Input: section "bulk" files like: LiveOcean_output/tef/[*]/bulk/[sect name].nc
where [*] = cas6_v3_lo8b_2017.01.01_2017.12.31 e.g.

Output: LiveOcean_output/tef/[*]/flux/two_layer_[season].p which is a pickled DataFrame whose index is the section names, and whose columns are: ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'].  Here _s and _f indicate that the layer is salty or fresh.  Also we have organized all the fluxes to be positive Eastward or Northward.  The averaging is over three "seasons" being:
- full=annual
- winter=JFM
- spring=AMJ
- summer=JAS
- fall=OND
and these are set by flux_fun.get_dtr(year)

------------------------------------------------------------------

(*) flux_plot_sections.py makes a plot of transport and salinity at all sections, along with a map, using the results of flux_make_two_layer.py above.

Input: LiveOcean_output/tef/[*]/flux/two_layer_[season].p

Output: LiveOcean_output/tef/tef_all_sections/all_sections_[year]_[season].png

..................................

(*) flux_plot_sections_clean.py is  a version of this that makes plots for the first paper.

Output: LiveOcean_output/tef/tef_all_sections_clean/all_sections_[year]_[season].png

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