README for x_moor

* mooring_extractor.py is the main driver for mooring extractions.  It can get a single mooring specified from the command line, or a list of moorings (and custom sets of variables) specified in moor_lists.py.

Input: LiveOcean ROMS history files

Output: A NetCDF file of the extracted mooring, typically packed as varname(time, z).

======================================================================
* moor_lists.py: Where I specify station names and locations as dict entries.  This is helpful when doing jobs with many moorings at once.

======================================================================
* user_moor_lists.py: This is not part of the repo, but instead is a placeholder for other users to make their own lists of moorings to extract.  Just copy moor_list.py to user_moor_list.py and edit it to define your own job_name.

======================================================================
* plot_mooring.py: makes a generic plot of a user-selected mooring extraction.

======================================================================
* add_C_vars.py: will add PH and ARAG to a user-selected mooring extraction.