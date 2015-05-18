function values = nc_getlast(ncfile, var, num_datums)
% NC_GETLAST:  Retrieves records at the end of an unlimited netCDF file
%
% DATA = NC_GETLAST(NCFILE,VARNAME,NUM_DATUMS) retrieves NUM_DATUMS 
% datums from the netCDF variable VARNAME in the netCDF file NCFILE.
% If NUM_DATUMS is not supplied, the default value is 1.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_getlast.m 2992 2010-03-30 17:22:26Z johnevans007 $
% $LastChangedDate: 2010-03-30 13:22:26 -0400 (Tue, 30 Mar 2010) $
% $LastChangedRevision: 2992 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~ischar(ncfile) 
	error ( 'SNCTOOLS:NC_GETLAST:badInput', ...
        'The netCDF file argument must be char.' );
end

if ~ischar(var) 
	error ( 'SNCTOOLS:NC_GETLAST:badInput', ...
        'The netCDF variable argument must be char.' );
end

if ( nargin == 2 )
	num_datums = 1;
else
	if ~isnumeric(num_datums) 
	    error ( 'SNCTOOLS:NC_GETLAST:badInput', ...
            'The num_datums argument must be numeric.' );
	end
	if num_datums <= 0
	    error ( 'SNCTOOLS:NC_GETLAST:badInput', ...
            'The num_datums argument must be positive.' );
	end

end

varlist = { var };
nb = nc_getbuffer ( ncfile, varlist, -1, num_datums );

values = nb.(var);

return

