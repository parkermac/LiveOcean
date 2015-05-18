function [S] = Z_get_S(infile)
% 1/30/2013  Parker MacCready
% written to use SNCTOOLS
%
% 'infile' is a ROMS history file
%
% this returns the structure:
%   S has vertical S-coordinate information
%
% added S.Vtransform 12/19/2012

% get vertical sigma-coordinate information (column vectors, bottom to top)
S.s_rho = nc_varget(infile,'s_rho');		% s-coordinate for "rho" (box mid-points)
S.s_w = nc_varget(infile,'s_w');         % s-coordinate for "w" (box interfaces)
S.hc = nc_varget(infile,'hc');			% s-coordinate critical depth (m)
S.Cs_r = nc_varget(infile,'Cs_r');       % s=coordinate structure function (rho)
S.Cs_w = nc_varget(infile,'Cs_w');       % s=coordinate structure function (w)
S.N = length(S.s_rho);        % number of vertical s-levels
S.Vtransform = nc_varget(infile,'Vtransform'); % Vtransform (e.g. 1 or 2)


