README for the flux_ code

This code is focused on the "efflux_reflux" concepts of Cokelet and Stewart (1985, and subsequent).  It picks of the TEF analysis from the point where we have done all the TEF sections, averaged into "bulk" inflows and outflows, and then forced them to be TWO-LAYER flows averaged over some time period like a year.  At that point we start working toward a super-simple numerical model in which the tracer concentation in "segments" (the volume in between two or more "sections").

------------------------------------------------------------------

* flux_fun.py is a module used throughout this code.

A key piece is the dict flux_fun.segs where we set up the segment names and define (i) the sections on each direction NSEW, and (ii) the rivers flowing into each segment.

Also has channel_dict and seg_dict which associates lists of sections or segments with each of four "chanels"

flux_fun.make_dist(x,y) to make vectors of distance along lon,lat vectors x,y

------------------------------------------------------------------

* flux_make_two_layer.py creates the time-mean two layer versions of the TEF variables, for each section.

Input: section "bulk" files like:
 - LiveOcean_output/tef/cas6_v3_lo8b/bulk/ai1.p

Output:

------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
------------------------------------------------------------------
