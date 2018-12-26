These programs create the ocean files used to force a roms run, typically for a single day or a 3 day forecast.  They bring in hycom data, either archived extractions for backfill, or from a web source for a forecast.  This is a fairly complicated process, becasue of some lack of regularity of the input formats.

===============================================================================
* make_forcing_main.py is the main driver, similar in basic construction and usage to all the code of the same name in other forcing folders (e.g. riv2, tide1, atm).

Input: command line arguments tell it what to do - specifically which gridname_tag = [gtag] to work on, what day [f_string], and whether this is a forecast or backfill.

Output: LiveOcean_output/[gtag]/[f_string]/ocn3/ocean_[clm, ini, bry].nc