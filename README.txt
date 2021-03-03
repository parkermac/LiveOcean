Parker MacCready - UW School of Oceanography

README for all the LiveOcean code.  Many of the sub-folders have their own detailed README's as well.  The overall purpose of this code is to run the LiveOcean daily forecast ROMS model.  It also supports any other ROMS model, such as our idealized estuary, with a few restrictions.

The code handles preprocessing to create the forcing files for a ROMS run, the drivers to automate making a run, and numerous post-processing tasks for after a run has finished.

Some assumptions that are baked into the code:

+ We are running ROMS.

+ Grids are "plaid", generated from vectors of longitude and latitude with np.meshgrid(lon, lat).  Stretching along either direction is allowed to increase resolution locally.

+ Runs are done in 1-day chunks and stored in day folders with names like "f2020.05.01" (YYYY.MM.DD), and so far we always store hourly output.  A given day folder may, however, have up to 73 hours of files if it is a three-day forecast.

+ Runs are typically named with a three-part word "[gridname]_[tag]_[ex_name]" such as cas6_v3_lo8b.  These codify the grid, the version of the forcing, and the version of the ROMS executable.

+ The software architecture is designed to be as modular as possible.  For example, to do a run with a different version of ROMS, you would only need to change the ex_name.

+ In general, all the drivers and primary code accept fairly consistent command line arguments.  This allows the code to be run for different cases and on different machines without editing.
An example of running the ocn4 driver shell-script on boiler:
  ./driver_forcing2.sh -g cas6 -t v3 -f ocn4 -r forecast > ./dlog_638_ocn4_a &
An example of running the plotting code in ipython:
  run pan_plot.py -g sj0 -t v0 -x lo8nest -0 2019.08.03
  
Associated repos, also on GitHub in parkermac:
* ptools (especially pgrid, the grid generation software)
* LiveOcean_roms/makefiles - centralizes and organizes the code for making ROMS executables
* LiveOcean_roms/LO_ROMS (private repo) - the version of ROMS with our NPZDO+Carbon code

Other required folders that are not repos:
* LiveOcean_data - has things like a coastline file, river climatologies, grids, HYCOM extractions, etc.  Typically these are large binary files that don't change often.

** Wherever this code exists it is assumed that these directories exist inside the parent directory:
LiveOcean (this repo)
LiveOcean_data
LiveOcean_output (where forcing files and post-processing products end up)
LiveOcean_roms

** To help with other people using this code, especially for post-processing tasks like particle tracking and extractions, I will place stubs in strategic places so that the code will first look for the corresponding file in a directory called "LiveOcean_user" (at the same level as LiveOcean) with corresponding sub-folders as needed.

** At a minimum a user needs to:
(1) Clone the LiveOcean repo to their machine(s), and do "git pull" periodically.
(2) Create LiveOcean_user/alpha/get_lo_info.sh with their own paths defined for each machine they will use.  LiveOcean_user is expected to be at the same level as LiveOcean in any directory structure.  They should make LiveOcean_user a repo they own in GitHub, and clone it as needed.
(3) Then for specific tasks they can create more files, mirroring the LiveOcean directory structure.  Presently only a few pieces of code are wired to look for versions in LiveOcean_user.  Please email me if you want a new connection point p.maccready@gmail.com.

-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------

* alpha/ contains top level code used by many programs.  Of greatest importance are get_lo_info.sh and Lfun.py.  There are also some of my own utility modules:
- zfun is for common tasks like tidal filtering
- zrfun is full of functions specific to working with ROMS files

=====================================================================================

* driver/ has all the drivers, mainly shell scripts (.sh) that you use to control pre- and post-processing operations.

=====================================================================================

* forcing/ has sub-folders for each of the individual pre- and post-processing jobs.  These are quite modular, and can be introduced to the forcing for a ROMS run just by changing a line in the dot_in file.  Each type of forcing (e.g. ocn4) typically has its own README and the same architecture: a primary program called make_forcing_main.py, and sometimes a subsidiary program, make_forcing_worker.m if MATLAB is needed.

The folder forcing/dot_in has sub-folders, one for each run (e.g. cas6_v3_lo8b) where the code and template for generating the ROMS .in file for a given day are kept.

Most of the code in forcing shares some input ant output structure that is contained in forcing_functions.py.  For example, this is the one place where we define acceptable command line inputs for each of the make_forcing_main.py. programs.

=====================================================================================

* plotting/ has code for plotting ROMS output.  The main program you use is pan_plot.py which accepts many arguments to decide which run to plot, what plotting code to use, time range, movie output, etc.  Needs a README!

=====================================================================================

* shared/ has folders of MATLAB modules like the seawater routines.

=====================================================================================

* superplot/ has code for the specific plotting task of making the year-long public-friendly movies on the LO website.

=====================================================================================

* tracker/ is the current particle tracking code, which is significantly faster becasue of extensive use of nearest-neighbor search using pre-created search trees (was called tracker2 until 2021.03.03).

=====================================================================================

* the x_[*]/ folders are all specific types of extraction code, meant to operate on existing ROMS output, or the piles of forcing files that went into making them.  They write their output to LiveOcean_output/[*].

- x_cast/ can very quickly make CTD cast-like extractions at a list of times and places.  A bit kludgy because it is hard-wired to look in some of my ptools_data directories to find observational records in a specific format that it then uses for locations and dates to extract.

- x_layer/ is a nice generic tool for quickly extracting layers over many history files and putting them in one NetCDF file.

- x_misc/ is not-so-generic extraction tools I keep for specific jobs.

- x_moor/ is for mooring extractions.  Has good capability to accept user-defined lists.

- x_river/ extracts river time series from LO forcing files

- x_svar/ is for salinity variance-related extractions.

- x_tef/ is for TEF extractions and subsequent analysis.  Currently (2020.05.02) very well organized but has a lot of ongoing development because it is central to a paper I am writing.

=====================================================================================
=====================================================================================
=====================================================================================
=====================================================================================
=====================================================================================

