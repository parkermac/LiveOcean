function [retrieval_method,fmt] = snc_read_backend(ncfile)
%SNC_READ_BACKEND   determine which netCDF library to use
%
% which backend do we employ?  Many, many possibilities to consider here.
%
%   [retrieval_method,fmt] = snc_read_backend(ncfile)
%
% returns selection for specified file, http or url.
%
%   [retrieval_method,fmt] = snc_read_backend()
%
% returns all available retrieval_method and fmt options
%
%See also: snctools, snc_format

retrieval_methods.java     = 'java';
retrieval_methods.tmw_hdf4 = 'tmw_hdf4';
retrieval_methods.mexnc    = 'mexnc';
retrieval_methods.tmw      = 'tmw';

fmts = snc_format();
fmt  = '';

if nargin==0
   retrieval_method = retrieval_methods;
else   
   retrieval_method = '';
end

use_java  = getpref('SNCTOOLS','USE_JAVA' ,false);
use_mexnc = getpref('SNCTOOLS','USE_MEXNC',false);

% Check for this early.
if isa(ncfile,'ucar.nc2.NetcdfFile') && use_java
    retrieval_method = retrieval_methods.java;
	fmt = fmts.netcdf_java;
	return
end


fmt = snc_format(ncfile);

% These cases have no alternatives.
if strcmp(fmt,fmts.HDF4) 
    % always use MATLAB's HDF interface for HDF-4 files.
    retrieval_method = retrieval_methods.tmw_hdf4;
	return
elseif use_java && (strcmp(fmt,fmts.GRIB) || strcmp(fmt,fmts.GRIB2) || strcmp(fmt,fmts.URL))
    % Always use netcdf-java for grib files or URLs (when java is enabled).
    retrieval_method = retrieval_methods.java;
	return
elseif strcmp(fmt,fmts.URL)
    % If java is not available, we have to assume that mexnc was compiled with
	% opendap support.
    retrieval_method = retrieval_methods.mexnc;
	return
end

mv = version('-release');
switch ( mv )
    case { '11', '12', '13' };
		error('Not supported on releases below R14.');

    case { '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
		% No native matlab support here.  We will favor java over
		% mexnc for now.
        if use_java && (strcmp(fmt,fmts.NetCDF) || strcmp(fmt,fmts.NetCDF4))
			% We will favor java over mexnc.
            retrieval_method = retrieval_methods.java;
        elseif use_mexnc && strcmp(fmt,fmts.NetCDF) 
            retrieval_method = retrieval_methods.mexnc;
        elseif use_mexnc && strcmp(fmt,fmts.NetCDF4) && use_mexnc
			% Assume the user knows what they are doing here.
        elseif use_java
            % Last chance is if it is some format that netcdf-java can handle.
            retrieval_method = retrieval_methods.java;
            fmt = fmts.netcdf_java;
        end
        
    case { '2008b', '2009a', '2009b', '2010a' }
        % 2008b introduced native netcdf-3 support.
		% netcdf-4 still requires either mexnc or java, and we will favor
		% java again.
        if strcmp(fmt,fmts.NetCDF) && use_mexnc
            % Use mexnc only if the user is dead serious about it.
            retrieval_method = retrieval_methods.mexnc;
        elseif strcmp(fmt,fmts.NetCDF)
            % otherwise use TMW for all local NC3 files 
            retrieval_method = retrieval_methods.tmw;
        elseif strcmp(fmt,fmts.NetCDF4) && use_java
            % Use TMW for all local NC3 files 
            retrieval_method = retrieval_methods.java;
        elseif strcmp(fmt,fmts.NetCDF4) 
            % Assume the user knows what they are doing.
            retrieval_method = retrieval_methods.mexnc;
        elseif use_java
            % Last chance is if it is some format that netcdf-java can handle.
            retrieval_method = retrieval_methods.java;
            fmt = fmts.netcdf_java;
        end

    otherwise
        % R2010b:  introduced netcdf-4 support.
        if strcmp(fmt,fmts.NetCDF) || strcmp(fmt,fmts.NetCDF4)
            retrieval_method = retrieval_methods.tmw;
        elseif use_java
            % Last chance is if it is some format that netcdf-java can handle.
            retrieval_method = retrieval_methods.java;
            fmt = fmts.netcdf_java;
        end

end

if isempty(retrieval_method)
    error('SNCTOOLS:unknownBackendSituation', ...
      'Could not determine which backend to use with %s.', ...
       ncfile );
end
return