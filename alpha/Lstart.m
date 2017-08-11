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
%
% it relies on the existence of a text file RUN_INFO.csv
% in which each line has two strings: an item name and its value
% separated by a comma

Ldir.gridname = gridname;
Ldir.tag = tag;
Ldir.gtag = [Ldir.gridname,'_',Ldir.tag];

% and get the parent
which_home = getenv('HOME');
switch which_home
    case '/Users/pm7';
        Ldir.env = 'pm_mac'
        Ldir.parent = [which_home,'/Documents/'];
    case '/home/parker'
        Ldir.env = 'fjord';
        Ldir.parent = '/data1/parker/';
    otherwise
        disp('Trouble filling out environment variables in Ldir')
end

%% set locations of things
Ldir.home = [Ldir.parent,'LiveOcean/'];
Ldir.out = [Ldir.parent,'LiveOcean_output/'];
Ldir.data = [Ldir.parent,'LiveOcean_data/'];

% Paths to shared code assumed to be available by many programs
addpath([Ldir.home,'shared/mexcdf/mexnc']);
addpath([Ldir.home,'shared/mexcdf/snctools']);
addpath([Ldir.home,'shared/seawater']);
addpath([Ldir.home,'shared/Z_functions']);


