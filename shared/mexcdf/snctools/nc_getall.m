function data = nc_getall ( ncfile )
%NC_GETALL Read entire contents of netCDF file.
%   This function is deprecated and is no longer supported.
% 
%   NCDATA = NC_GETALL(NCFILE) reads the entire contents of the netCDF 
%   file NCFILE into the structure NCDATA.  


warning ( 'SNCTOOLS:nc_getall:deprecated', ...
          'NC_GETALL is deprecated and may be removed in a future version of SNCTOOLS.');


% Show usage if too few arguments.
if nargin~=1 
    error ( 'must have one input argument.\n' );
end


switch ( version('-release') )
    case { '11', '12', '13', '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
        data = nc_getall_mex ( ncfile );
    otherwise
        data = nc_getall_tmw ( ncfile );
end


%-----------------------------------------------------------------------
function data = nc_getall_mex ( ncfile )
data = [];

[cdfid,status ]=mexnc('open',ncfile,'NOWRITE');
if status ~= 0
    error ( mexnc('strerror', status) );
end

[nvars, status] = mexnc('INQ_NVARS', cdfid);
if status < 0
    mexnc('close',cdfid);
    error ( mexnc('strerror', status) );
end
[ngatts, status] = mexnc('INQ_NATTS', cdfid);
if status < 0
    mexnc('close',cdfid);
    error ( mexnc('strerror', status) );
end

for varid=0:nvars-1

    varstruct = [];

    [varname, datatype, ndims, dims, natts, status] = mexnc('INQ_VAR', cdfid, varid); %#ok<ASGLU>
    if status < 0 
        mexnc('close',cdfid);
        error ( mexnc('strerror', status) );
    end

    % If ndims is zero, then it must be a singleton variable.  Don't bother trying
    % to retrieve the data, there is none.
    if ( ndims == 0 )
        varstruct.data = [];
    else
        values = nc_varget(ncfile, varname);
        varstruct.data = values;
    end

    % get all the attributes
    for attnum = 0:natts-1

        [attname, status] = mexnc('inq_attname', cdfid, varid, attnum);
        if status < 0 
            mexnc('close',cdfid);
            error ( mexnc('strerror', status) );
        end

        try
            attval = nc_attget(ncfile, varname, attname);
        catch %#ok<CTCH>
            mexnc('close',cdfid);
            error ( 'SNCTOOLS:nc_getall:mexnc:attributeRetrievalFailed', ...
                'Failed to retrieve attribute #%d, ''%s''', ...
                attnum, lasterr); %#ok<LERR>
        end
        
        % Matlab structures don't like the leading '_'
        if strcmp(attname,'_FillValue' )
            attname = 'FillValue';
        end

        sanitized_attname = genvarname(attname);

        % this puts the attribute into the variable structure
        varstruct.(sanitized_attname) = attval;

    end


    % Add this variable to the entire file structure
    data.(varname) = varstruct;

end


% Do the same for the global attributes
global_atts = [];
for attnum = 0:ngatts-1

    [attname, status] = mexnc('inq_attname', cdfid, nc_global, attnum);
    if status < 0 
        mexnc('close',cdfid);
        error ( mexnc('strerror',status) );
    end

    try
        attval = nc_attget(ncfile, nc_global, attname);
    catch %#ok<CTCH>
        mexnc('close',cdfid);
        error ( 'SNCTOOLS:nc_getall:mexnc:globalAttRetrievalFailed', ...
            'Failed to retrieve global attribute #%d, ''%s''', ...
            attnum, lasterr); %#ok<LERR>
    end
    
    sanitized_attname = genvarname(attname);

    % this puts the attribute into the variable structure
    global_atts.(sanitized_attname) = attval;

end

if ~isempty ( global_atts )
    data.global_atts = global_atts;
end

mexnc('close',cdfid);

if isempty(data)
    data = struct([]);
end

return




