function [G,S,T] = Z_get_basic_info(infile)
% 12/19/2012  Parker MacCready
% written to use SNCTOOLS
%
% 'infile' is a ROMS history file
%
% this returns the structures:
%   G has horizontal grid info including bathymetry
%   S has vertical S-coordinate information
%   T has time information
%
% fields are packed as (t,k,j,i) with k increasing upwards
% and we are assuming that there is only one time level per file
%
% added S.Vtransform 12/19/2012

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

% get vertical sigma-coordinate information (column vectors, bottom to top)
S.s_rho = nc_varget(infile,'s_rho');		% s-coordinate for "rho" (box mid-points)
S.s_w = nc_varget(infile,'s_w');         % s-coordinate for "w" (box interfaces)
S.hc = nc_varget(infile,'hc');			% s-coordinate critical depth (m)
S.Cs_r = nc_varget(infile,'Cs_r');       % s=coordinate structure function (rho)
S.Cs_w = nc_varget(infile,'Cs_w');       % s=coordinate structure function (w)
S.N = length(S.s_rho);        % number of vertical s-levels
S.Vtransform = nc_varget(infile,'Vtransform'); % Vtransform (e.g. 1 or 2)

% time info
T.ocean_time = nc_varget(infile,'ocean_time'); % time in seconds
T.dstart = nc_varget(infile,'dstart'); % day of the first save of this run
T.dstart_units = nc_attget(infile,'dstart','units');
% NOTE T.dstart_units is a string giving the time reference
% e.g. seconds since 2004-01-01 00:00:00
iii = find(T.dstart_units == '-');
iii = iii(1);
year0 = str2num(T.dstart_units(iii-4:iii-1));
month0 = str2num(T.dstart_units(iii+1:iii+2));
day0 = str2num(T.dstart_units(iii+4:iii+5));
iii = find(T.dstart_units == ':');
iii = iii(1);
hour0 = str2num(T.dstart_units(iii-2:iii-1));
minute0 = str2num(T.dstart_units(iii+1:iii+2));
second0 = str2num(T.dstart_units(iii+4:iii+5));
datenum0 = datenum(year0,month0,day0,hour0,minute0,second0);
T.time_datenum = datenum0 + T.ocean_time/86400;


