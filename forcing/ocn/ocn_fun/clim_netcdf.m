function [V,AAlon,AAlat,AAmask] = ...
    clim_netcdf(indir,varname,gridfile)
% 4/27/2012  Parker MacCready
%
% returns variable-specific information

% Read in the grid matrices
lon_rho = nc_varget(gridfile,'lon_rho');
lat_rho = nc_varget(gridfile,'lat_rho');
mask_rho = nc_varget(gridfile,'mask_rho');

lon_u = nc_varget(gridfile,'lon_u');
lat_u = nc_varget(gridfile,'lat_u');
mask_u = nc_varget(gridfile,'mask_u');

lon_v = nc_varget(gridfile,'lon_v');
lat_v = nc_varget(gridfile,'lat_v');
mask_v = nc_varget(gridfile,'mask_v');

do_depth = 1; % flat to indicate that there is depth
do_addtimevar = 1; % include the time variable (set to 0 for v and vbar)

switch varname
    case 's3d'
        invarname = 's3d_filt';
        infile = [indir,varname,'.nc'];
        ncvarname = 'salt';
        nclongname = 'salinity climatology';
        ncunits = 'PSU';
        nctimename = 'salt_time';
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
    case 't3d'
        invarname = 't3d_filt';
        infile = [indir,varname,'.nc'];
        ncvarname = 'temp';
        nclongname = 'potential temperature climatology';
        ncunits = 'Celsius';
        nctimename = 'temp_time';
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
    case 'u3d'
        invarname = 'u3d_filt';
        infile = [indir,varname,'.nc'];
        ncvarname = 'u';
        nclongname = 'u-momentum component climatology';
        ncunits = 'meter second-1';
        nctimename = 'v3d_time';
        ncxiname = 'xi_u'; ncetaname = 'eta_u';
        AAlon = lon_u; AAlat = lat_u; AAmask = mask_u;
    case 'v3d'
        invarname = 'v3d_filt';
        infile = [indir,varname,'.nc'];
        ncvarname = 'v';
        nclongname = 'v-momentum component climatology';
        ncunits = 'meter second-1';
        nctimename = 'v3d_time';
        do_addtimevar = 0;
        ncxiname = 'xi_v'; ncetaname = 'eta_v';
        AAlon = lon_v; AAlat = lat_v; AAmask = mask_v;
    case 'ssh'
        invarname = 'ssh_filt';
        infile = [indir,varname,'.nc'];
        ncvarname = 'zeta';
        nclongname = 'sea surface height climatology';
        ncunits = 'meter';
        nctimename = 'zeta_time';
        do_depth = 0; % 2D variable
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
    case 'ubar'
        invarname = 'ubar';
        infile = [];
        ncvarname = 'ubar';
        nclongname = 'vertically averaged u-momentum climatology';
        ncunits = 'meter second-1';
        nctimename = 'v2d_time';
        do_depth = 0; % 2D variable
        ncxiname = 'xi_u'; ncetaname = 'eta_u';
        AAlon = lon_u; AAlat = lat_u; AAmask = mask_u;
    case 'vbar'
        invarname = 'vbar';
        infile = [];
        ncvarname = 'vbar';
        nclongname = 'vertically averaged v-momentum climatology';
        ncunits = 'meter second-1';
        nctimename = 'v2d_time';
        do_depth = 0; % 2D variable
        do_addtimevar = 0;
        ncxiname = 'xi_v'; ncetaname = 'eta_v';
        AAlon = lon_v; AAlat = lat_v; AAmask = mask_v;
end

% pack the netcdf info in a structure "V"
V.invarname = invarname; V.infile = infile; V.ncvarname = ncvarname;
V.nclongname = nclongname; V.ncunits = ncunits; V.nctimename = nctimename;
V.do_depth = do_depth; V.do_addtimevar = do_addtimevar;
V.ncxiname = ncxiname; V.ncetaname = ncetaname;
