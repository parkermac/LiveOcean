function [] = clim_netcdf_new(gridfile,outfile,N)
% 12/6/2011  Parker MacCready
%
% creates the initial climatology file

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

[Mrho,Lrho] = size(lon_rho);
[Mu,Lu] = size(lon_u);
[Mv,Lv] = size(lon_v);

% create the NetCDF output file
my_mode = bitor (nc_clobber_mode, nc_64bit_offset_mode );
nc_create_empty(outfile, my_mode );
% global attributes
nc_padheader (outfile, 20000 );
nc_attput(outfile, nc_global, 'type','ROMS Climatology File');
% define some dimensions
nc_add_dimension(outfile, 'xi_rho', Lrho);
nc_add_dimension(outfile, 'eta_rho', Mrho);
nc_add_dimension(outfile, 'xi_u', Lu);
nc_add_dimension(outfile, 'eta_u', Mu);
nc_add_dimension(outfile, 'xi_v', Lv);
nc_add_dimension(outfile, 'eta_v', Mv);
nc_add_dimension(outfile, 's_rho', N);

