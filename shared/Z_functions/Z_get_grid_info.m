function [G] = Z_get_grid_info(infile)
% 3/2/2011  Parker MacCready
% written to use SNCTOOLS
%
% 'infile' is a ROMS history file
%
% this returns the structures:
%   G has horizontal grid info including bathymetry
%

% get grid and bathymetry
G.h = nc_varget(infile,'h');             % depth (m) of bathymetry (positive down)
G.lon_rho = nc_varget(infile,'lon_rho'); G.lat_rho = nc_varget(infile,'lat_rho');
G.lon_u = nc_varget(infile,'lon_u'); G.lat_u = nc_varget(infile,'lat_u');
G.lon_v = nc_varget(infile,'lon_v'); G.lat_v = nc_varget(infile,'lat_v');
G.lon_psi = nc_varget(infile,'lon_psi'); G.lat_psi = nc_varget(infile,'lat_psi');
G.mask_rho = logical(nc_varget(infile,'mask_rho'));
G.mask_u = nc_varget(infile,'mask_u'); G.mask_v = nc_varget(infile,'mask_v');
pm = nc_varget(infile,'pm'); pn = nc_varget(infile,'pn');
G.DX = 1./pm; G.DY = 1./pn; % grid sizes (m)
[G.M,G.L] = size(G.h);
