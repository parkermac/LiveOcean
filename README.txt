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
An example of running the ocn4 driver shellscript on boiler:
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

* alpha/ contains top level code used by many programs.  In...

=====================================================================================
=====================================================================================
=====================================================================================
=====================================================================================

* the x_[*]/ folders are all specific extraction, meant to operate on existing ROMS output, or the piles of forcing files that went into making them.  They write their output to LiveOcean_output/[*].

=====================================================================================
=====================================================================================
=====================================================================================
=====================================================================================
=====================================================================================

