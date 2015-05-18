function [coast_handle] = Z_addcoast(coastfile)
% 1/25/2014  Parker MacCready
%
% Adds a coastline to an existing plot.  Without an explicit coastfile
% it relies on the path to default being set in Lstart.m.
%
% Usage:
%
% [coast_handle] = Z_addcoast(coastfile);
%
% or just:
%
% Z_addcoast;

if nargin == 0
    coastfile = 'pnw_coast_combined.mat';
end
hold on
load(coastfile);
coast_handle = plot(lon,lat,'-k');
hold off

