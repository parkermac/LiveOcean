README for x_moor

=====================================================================
* mooring_extractor.py is the main driver for mooring extractions.  It can get a single mooring specified from the command line, or a list of moorings (and custom sets of variables) specified in moor_lists.py or LiveOcean_user/x_moor/user_moor_lists.py.

Input: LiveOcean ROMS history files

Output: A NetCDF file of the extracted mooring, typically packed as varname(time, z).

NOTE: mooring_extractor_fast.py is development code - not likely to be useful.

======================================================================
* moor_lists.py: Where I specify station names and locations as dict entries.  This is helpful when doing jobs with many moorings at once.

======================================================================
* LiveOcean_user/x_moor/user_moor_lists.py: This is not part of the repo, but instead is a placeholder for other users to make their own lists of moorings to extract.  Just copy moor_lists.py to user_moor_lists.py and edit it to define your own job_name(s).

======================================================================
* plot_mooring.py: makes a generic plot of a user-selected mooring extraction.

======================================================================
* add_C_vars.py: will add PH and ARAG to a user-selected mooring extraction.

======================================================================
* add_C_vars_multifile.py: will add PH and ARAG to  all the mooring extractions
in a user-selected directory.