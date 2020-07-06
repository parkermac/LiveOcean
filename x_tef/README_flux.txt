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

Output: LiveOcean_output/tef/[*]/flux/two_layer_[season].p which is a pickled DataFrame whose index is the section names, and whose columns are: ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'].  Here _s and _f indicate that the layer is salty or fresh.  Also we have organized all the fluxes to be positive Eastward or Northward.  The averaging is over three "seasons" being:
- full=annual
- winter=JFM
- spring=AMJ
- summer=JAS
- fall=OND
and these are set by flux_fun.get_dtr(year)

------------------------------------------------------------------

* flux_plot_sections.py makes a plot of transport and salinity at all sections, along with a map, using the results of flux_make_two_layer.py above.

Input: LiveOcean_output/tef/[*]/flux/two_layer_[season].p

Output: LiveOcean_output/tef/tef_all_sections/all_sections_[year]_[season].png

..................................

* flux_plot_sections_clean.py is  a version of this that makes plots for the first paper.

Output: LiveOcean_output/tef/tef_all_sections_clean/all_sections_[year]_[season].png

------------------------------------------------------------------

* flux_get_vol.py gets the volume (with SSH = 0) of each segment.  Uses a cute "mine sweeper" like algorithm to find all the unmasked rho-grid cells on the grid between the sections that define a segment.

Input: flux_fun.segs, the model grid (hard wired in the code to be cas6), and sect_df = tef_fun.get_sect_df().

Output: LiveOcean_output/tef/volumes_cas6/...
	volumes.p a pickled DataFrame with segment names for the index and columns: ['volume m3', 'area m2', 'lon', 'lat']
	bathy_dict.p pickled dict of the grid bathy, for plotting
	ji_dict.p pickled dict (keys = segment names) of the ji indices of each segment area, also for plotting
	
* flux_V_plot.py makes a simple plot of the segment volumes represented as rectangles:
	LiveOcean_output/tef/misc_figs_cas6/volume_plot.png
	Eventually this might be an interesting framework for plotting the transports used in the flux engine.
	
------------------------------------------------------------------

* flux_seg_sect_maps.py

Input: grid, rivers, sections and segments (from volumes.py)

Output: LiveOcean_output/tef/misc_figs_cas6/seg_map_full.png and seg_map_focus.py which are useful for seeing how each segment is related to its defining sections and rivers.  Must be consistent with flux_fun.segs.

------------------------------------------------------------------

* flux_seg_map_simple.py

Input: the output of flux_get_vol.py

Output: LiveOcean_output/tef/misc_figs_cas6/seg_map_simple.png, a graphically appealing plot showing the four basins defined by our "channels" as a colorful map.

------------------------------------------------------------------

* flux_get_A.py is a central and complex part of the flux machinery.  Its job is to take the TEF transports created by flux_make_two_layer.py, and the river flows and turn them into efflux-reflux fractions used by the flux engine.

Input: LiveOcean_output/tef/[*]/flux/two_layer_[season].p and...
	river flow extracted over the model year from the forcing (see x_river)
	flux_fun.segs

Output: LiveOcean_output/tef/[*]/flux/q_df_[season].p a pickled DataFrame where the index is the segment list, but with separate salty and fresh entries, like J1_s and J1_f.  The columns are the same, but also there are four more columns at the start: 'ocean_s', 'ocean_f', 'river_s', 'river_f'. From the comments:
# The way q_df is used in the flux_engine is that the row denotes where net transports are
# ending up, and the columns denote where the transports are coming from.  Each entry in
# q_df is a transport (m3/s).  Hence when we multiply a row by the tracer concentrations
# in each segment and then sum along that row, we get the net flux of tracer into or out of
# the segment specified by that row.

NOTE: Look deep in the code for methods of calculation, and hand adjustments to, the efflux-reflux coefficients.

------------------------------------------------------------------

* flux_vtrans_plot.py makes a plot of vertical transport in all the segments.

Input: volumes.p and q_df_[season].p created above

Output:LiveOcean_output/tef/vertical_transport_plots/vtrans_[year].png

------------------------------------------------------------------

* *** flux_engine.py ***
This is the main piece of code this whole efflux-reflux analysis has been driving toward.  It takes the transport matrix created by flux_get_A.py, and the volumes, and then uses them to do forward calculations of time-dependent tracer fields.  You can control the boundary conditions (where the tracer comes from), and whether it sinks - like organic particles, using command line arguments.

Input: volumes.p and q_df_[season].p created above

Output: indir = LiveOcean_output/tef/[*]/flux/
    cc.to_pickle(indir + 'cc_' + source + '_' + season + sink_tag + '.p') # FINAL STATE WITH AGE
    aa.to_pickle(indir + 'aa_' + source + '_' + season + sink_tag + '.p') # TIME DEPENDENT STATE
	
- cc is a DataFrame, same index as q_df, and the columns 'c' and 'ca' are the final (after say 6 years) values of the tracer, and an aging tracer.

 - aa is a DataFrame of a full time-dependent array, which you could use to calculate residence times in experiments where you are setting up an initial condition.  For aa the index is time, and the columns are the q_df indices.

------------------------------------------------------------------

* flux_plot_validation.py makes a plot of the steady-state output of the flux engine for the case which is supposed to reproduce the ocean salinity.  This is the only way to test the flux engine and its efflux-reflux coefficients.  It compares the cc_ values to the two_layer_[season].p TEF values.

Input: cc_ocean_salt_[season].p (from flux_engine.py) and two_layer_[season].p

Output:LiveOcean_output/tef/validation_plots/validation_plot_[year]_full.png

NOTE: This is  only meaningful for the full annual average, because it finds an equilibrated solution.
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

Output: LiveOcean_output/tef/[*]/flux/two_layer_[season].p which is a pickled DataFrame whose index is the section names, and whose columns are: ['q_s', 'q_f', 'f_s', 'f_f', 's_s', 's_f', 'lon', 'lat'].  Here _s and _f indicate that the layer is salty or fresh.  Also we have organized all the fluxes to be positive Eastward or Northward.  The averaging is over three "seasons" being:
- full=annual
- winter=JFM
- spring=AMJ
- summer=JAS
- fall=OND
and these are set by flux_fun.get_dtr(year)

------------------------------------------------------------------

* flux_plot_sections.py makes a plot of transport and salinity at all sections, along with a map, using the results of flux_make_two_layer.py above.

Input: LiveOcean_output/tef/[*]/flux/two_layer_[season].p

Output: LiveOcean_output/tef/tef_all_sections/all_sections_[year]_[season].png

..................................

* flux_plot_sections_clean.py is  a version of this that makes plots for the first paper.

Output: LiveOcean_output/tef/tef_all_sections_clean/all_sections_[year]_[season].png

------------------------------------------------------------------

* flux_get_vol.py gets the volume (with SSH = 0) of each segment.  Uses a cute "mine sweeper" like algorithm to find all the unmasked rho-grid cells on the grid between the sections that define a segment.

Input: flux_fun.segs, the model grid (hard wired in the code to be cas6), and sect_df = tef_fun.get_sect_df().

Output: LiveOcean_output/tef/volumes_cas6/...
	volumes.p a pickled DataFrame with segment names for the index and columns: ['volume m3', 'area m2', 'lon', 'lat']
	bathy_dict.p pickled dict of the grid bathy, for plotting
	ji_dict.p pickled dict (keys = segment names) of the ji indices of each segment area, also for plotting
	
* flux_V_plot.py makes a simple plot of the segment volumes represented as rectangles:
	LiveOcean_output/tef/misc_figs_cas6/volume_plot.png
	Eventually this might be an interesting framework for plotting the transports used in the flux engine.
	
------------------------------------------------------------------

* flux_seg_sect_maps.py

Input: grid, rivers, sections and segments (from volumes.py)

Output: LiveOcean_output/tef/misc_figs_cas6/seg_map_full.png and seg_map_focus.py which are useful for seeing how each segment is related to its defining sections and rivers.  Must be consistent with flux_fun.segs.

------------------------------------------------------------------

* flux_seg_map_simple.py

Input: the output of flux_get_vol.py

Output: LiveOcean_output/tef/misc_figs_cas6/seg_map_simple.png, a graphically appealing plot showing the four basins defined by our "channels" as a colorful map.

------------------------------------------------------------------

* flux_get_A.py is a central and complex part of the flux machinery.  Its job is to take the TEF transports created by flux_make_two_layer.py, and the river flows and turn them into efflux-reflux fractions used by the flux engine.

Input: LiveOcean_output/tef/[*]/flux/two_layer_[season].p and...
	river flow extracted over the model year from the forcing (see x_river)
	flux_fun.segs

Output: LiveOcean_output/tef/[*]/flux/q_df_[season].p a pickled DataFrame where the index is the segment list, but with separate salty and fresh entries, like J1_s and J1_f.  The columns are the same, but also there are four more columns at the start: 'ocean_s', 'ocean_f', 'river_s', 'river_f'. From the comments:
# The way q_df is used in the flux_engine is that the row denotes where net transports are
# ending up, and the columns denote where the transports are coming from.  Each entry in
# q_df is a transport (m3/s).  Hence when we multiply a row by the tracer concentrations
# in each segment and then sum along that row, we get the net flux of tracer into or out of
# the segment specified by that row.

NOTE: Look deep in the code for methods of calculation, and hand adjustments to, the efflux-reflux coefficients.

------------------------------------------------------------------

* flux_vtrans_plot.py makes a plot of vertical transport in all the segments.

Input: volumes.p and q_df_[season].p created above

Output:LiveOcean_output/tef/vertical_transport_plots/vtrans_[year].png

------------------------------------------------------------------

* *** flux_engine.py ***
This is the main piece of code this whole efflux-reflux analysis has been driving toward.  It takes the transport matrix created by flux_get_A.py, and the volumes, and then uses them to do forward calculations of time-dependent tracer fields.  You can control the boundary conditions (where the tracer comes from), and whether it sinks - like organic particles, using command line arguments.

Input: volumes.p and q_df_[season].p created above

Output: indir = LiveOcean_output/tef/[*]/flux/
    cc.to_pickle(indir + 'cc_' + source + '_' + season + sink_tag + '.p') # FINAL STATE WITH AGE
    aa.to_pickle(indir + 'aa_' + source + '_' + season + sink_tag + '.p') # TIME DEPENDENT STATE
	
- cc is a DataFrame, same index as q_df, and the columns 'c' and 'ca' are the final (after say 6 years) values of the tracer, and an aging tracer.

 - aa is a DataFrame of a full time-dependent array, which you could use to calculate residence times in experiments where you are setting up an initial condition.  For aa the index is time, and the columns are the q_df indices.

------------------------------------------------------------------

* flux_plot_validation.py makes a plot of the steady-state output of the flux engine for the case which is supposed to reproduce the ocean salinity.  This is the only way to test the flux engine and its efflux-reflux coefficients.  It compares the cc_ values to the two_layer_[season].p TEF values.

Input: cc_ocean_salt_[season].p (from flux_engine.py) and two_layer_[season].p

Output:LiveOcean_output/tef/validation_plots/validation_plot_[year]_full.png

NOTE: This is  only meaningful for the full annual average, because it finds an equilibrated solution.

------------------------------------------------------------------

* flux_plot_cc.py makes a plot of tracer concentration and age from any instance of a cc_ file created by flux_engine.py.

------------------------------------------------------------------

* flux_residence time.py makes a movie out of times saved in any instance of an aa_ file created by flux_engine.py.

Input: all of the "IC_" series of flux_engine calculations.

Output:
- Beautiful plot of the residence time: LiveOcean_output/tef/misc_figs_cas6/all_residence_times.png
- DataFrame of all residence times: LiveOcean_output/tef/misc_figs_cas6/tres_df.p, with columns "tres" [residence time in days] and "tres0" [residence time in days without reflux (V/Qin)]

------------------------------------------------------------------

* flux_V_movie.py makes a movie from one of the S_ or IC_ runs, including information about residence time (with and without reflux).  The graphical format is the series of volume rectangles as in flux_V_plot.py.

Input: any of the "IC_" series of flux_engine runs, e.g. [run] = IC_HoodCanalInner_2018_fall, and LiveOcean_output/tef/misc_figs_cas6/tres_df.p

Output: LiveOcean_output/tef/movies/[gtagex]/[run]/movie.mp4


------------------------------------------------------------------

* flux_get_s.py creates time series of hourly salinity and volume in segments, over a user specified (like a year) period.

Input: ROMS history files

Output: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p, pickled DataFrames where the index is hourly datetime and the columns are the segments.

------------------------------------------------------------------

* flux_lowpass_s.py creates time series of daily salinity, volume, and net salt in segments.

Input: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p

Output: LiveOcean_output/tef/[*]/flux/daily_segment_[salinity,volume,net_salt].p

------------------------------------------------------------------

* flux_salt_budget.py makes a complete volume and salt budget for a user-specified set of segments.
These budgets are of the form:
	dSnet/dt = QSin + QSout, and
	dV/dt = Qin + Qout + Qr
They are mainly useful for knowing that the budgets add up in each of the basins and years to reasonable accuracy.  The values plotted are DAILY (tidally averaged).

Input: LiveOcean_output/tef/[*]/flux/daily_segment_[volume,net_salt].p

Output: A screen plot of the budget vs. time, and a saved png such as:
LiveOcean_output/tef/salt_budget_plots/salt_budget_2018_Salish_Sea.png

------------------------------------------------------------------

* flux_qe_salt_budget.py is like flux_salt_budget.py above, but it is more focused on dynamics controlling things at a sill.  Like flux_salt_budget.py it plots daily (tidally averaged) values over a year.
It plots a rearranged salt budget of the form:
	dSnet/dt = -Qr*Sbar + Qe*DS 
where Qe is like (Qin-Qout)/2, and Sbar is like (Sin + Sout)/2, and these are exact recombinations of the terms in flux_salt_budget.py.

------------------------------------------------------------------

* flux_sill_dyn.py uses much of the same information as in flux_qe_salt_budget.py above, but the terms are for the segments in the middle of a sill like ai2, and then we calculate some quantities that MAY control the exchange flow salt flux (Qe*DS) at that segment, such as dSbar/dx (calculated from ai1 and ai3 for example), as well as the net tidal dissipation between those sections (calculated as the difference of net tidal energy flux) which is a proxy for Kv.  These various possible controls are plotted as scatterplots.  Works with DAILY values for a year.

------------------------------------------------------------------

* flux_3year_salt.py is also like flux_salt_budget.py above, but it seeks to shed light on what controls the variability of mean salinity in a chosen basin.  It plots MONTHLY averages of mean salinity in the basin over THREE YEARS, along with other info like the TEF terms and a salt budget.

------------------------------------------------------------------
* flux_plot_sections.py makes maps of year-averaged tidal energy flux and net transport
at all sections

Input: section "bulk" files

Output: ...LiveOcean_output/tef/misc_figs_cas6/Tide_and_Qnet_[year].png
------------------------------------------------------------------
------------------------------------------------------------------

------------------------------------------------------------------

* flux_plot_cc.py makes a plot of tracer concentration and age from any instance of a cc_ file created by flux_engine.py.

------------------------------------------------------------------

* flux_residence time.py makes a movie out of times saved in any instance of an aa_ file created by flux_engine.py.

Input: all of the "IC_" series of flux_engine calculations.

Output:
- Beautiful plot of the residence time: LiveOcean_output/tef/misc_figs_cas6/all_residence_times.png
- DataFrame of all residence times: LiveOcean_output/tef/misc_figs_cas6/tres_df.p, with columns "tres" [residence time in days] and "tres0" [residence time in days without reflux (V/Qin)]

------------------------------------------------------------------

* flux_V_movie.py makes a movie from one of the S_ or IC_ runs, including information about residence time (with and without reflux).  The graphical format is the series of volume rectangles as in flux_V_plot.py.

Input: any of the "IC_" series of flux_engine runs, e.g. [run] = IC_HoodCanalInner_2018_fall, and LiveOcean_output/tef/misc_figs_cas6/tres_df.p

Output: LiveOcean_output/tef/movies/[gtagex]/[run]/movie.mp4


------------------------------------------------------------------

* flux_get_s.py creates time series of hourly salinity and volume in segments, over a user specified (like a year) period.

Input: ROMS history files

Output: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p, pickled DataFrames where the index is hourly datetime and the columns are the segments.

------------------------------------------------------------------

* flux_lowpass_s.py creates time series of daily salinity, volume, and net salt in segments.

Input: LiveOcean_output/tef/[*]/flux/hourly_segment_[salinity,volume].p

Output: LiveOcean_output/tef/[*]/flux/daily_segment_[salinity,volume,net_salt].p

------------------------------------------------------------------

* flux_salt_budget.py makes a complete volume and salt budget for a user-specified set of segments.
These budgets are of the form:
	dSnet/dt = QSin + QSout, and
	dV/dt = Qin + Qout + Qr
They are mainly useful for knowing that the budgets add up in each of the basins and years to reasonable accuracy.  The values plotted are DAILY (tidally averaged).

Input: LiveOcean_output/tef/[*]/flux/daily_segment_[volume,net_salt].p

Output: A screen plot of the budget vs. time, and a saved png such as:
LiveOcean_output/tef/salt_budget_plots/salt_budget_2018_Salish_Sea.png

------------------------------------------------------------------

* flux_qe_salt_budget.py is like flux_salt_budget.py above, but it is more focused on dynamics controlling things at a sill.  Like flux_salt_budget.py it plots daily (tidally averaged) values over a year.
It plots a rearranged salt budget of the form:
	dSnet/dt = -Qr*Sbar + Qe*DS 
where Qe is like (Qin-Qout)/2, and Sbar is like (Sin + Sout)/2, and these are exact recombinations of the terms in flux_salt_budget.py.

------------------------------------------------------------------

* flux_sill_dyn.py uses much of the same information as in flux_qe_salt_budget.py above, but the terms are for the segments in the middle of a sill like ai2, and then we calculate some quantities that MAY control the exchange flow salt flux (Qe*DS) at that segment, such as dSbar/dx (calculated from ai1 and ai3 for example), as well as the net tidal dissipation between those sections (calculated as the difference of net tidal energy flux) which is a proxy for Kv.  These various possible controls are plotted as scatterplots.  Works with DAILY values for a year.

------------------------------------------------------------------

* flux_3year_salt.py is also like flux_salt_budget.py above, but it seeks to shed light on what controls the variability of mean salinity in a chosen basin.  It plots MONTHLY averages of mean salinity in the basin over THREE YEARS, along with other info like the TEF terms and a salt budget.

------------------------------------------------------------------
* flux_plot_sections.py makes maps of year-averaged tidal energy flux and net transport
at all sections

Input: section "bulk" files

Output: ...LiveOcean_output/tef/misc_figs_cas6/Tide_and_Qnet_[year].png
------------------------------------------------------------------
------------------------------------------------------------------
