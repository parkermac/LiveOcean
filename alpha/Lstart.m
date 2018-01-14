function [Ldir] = Lstart(gridname, tag)
% Parker MacCready
%
% A function to be invoked at the start of any primary MATLAB code anywhere
% in the LiveOcean folder.  It returns a structure Ldir that has handy
% pathnames, and it adds useful toolboxes to the path.
%
% Typical usage (depending on directory location in LiveOcean/):
%
% addpath('../alpha'); Ldir = Lstart(gridname, tag);

% find the path to alpha (assumes path has been added)
alp0 = which('Lstart');
alp = alp0(1:end-8);
% read lo_info into a struct
fid = fopen([alp,'lo_info.csv']);
while ~feof(fid)
    line = fgetl(fid);
    ic = find(line == ',');
    name = line(1:ic-1);
    value = line(ic+1:end);
    Ldir.(name) = value;
end
fclose(fid);

Ldir.gridname = gridname;
Ldir.tag = tag;
Ldir.gtag = [Ldir.gridname,'_',Ldir.tag];

%% set locations of things

% Paths to shared code assumed to be available by many programs
addpath([Ldir.LO,'shared/mexcdf/mexnc']);
addpath([Ldir.LO,'shared/mexcdf/snctools']);
addpath([Ldir.LO,'shared/seawater']);
addpath([Ldir.LO,'shared/Z_functions']);


