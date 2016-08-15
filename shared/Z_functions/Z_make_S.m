function Z_make_S(grid_dir)
%
% This makes and saves the structure "S" for a given grid.
% It is designed to be run as a subprocess of
% pgrid/grid_to_LiveOcean.py.  Mainly is it a wrapper
% for Z_scoord().

% This gives access to nc_varget.
addpath('../../alpha'); Ldir = Lstart('foo', 'foo');

outname = [grid_dir,'S.mat'];

% determine what hmin was used in the grid
gfile = [grid_dir,'grid.nc'];
h = nc_varget(gfile,'h');
hmin = max(0, min(h(:)));

% get the remaining S-coordinate parameters from a text file
fid = fopen([grid_dir,'S_COORDINATE_INFO.csv'],'r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
items = C{1};
values = C{2};
for ii = 1:length(items)
    if(strcmp(items{ii},'THETA_S'))
        theta_s = str2num(values{ii});
    elseif(strcmp(items{ii},'THETA_B'))
        theta_b = str2num(values{ii});
    elseif(strcmp(items{ii},'TCLINE'))
        tcline = str2num(values{ii});
    elseif(strcmp(items{ii},'N'))
        N = str2num(values{ii});
    elseif(strcmp(items{ii},'VTRANSFORM'))
        Vtransform = str2num(values{ii});
    elseif(strcmp(items{ii},'VSTRETCHING'))
        Vstretching = str2num(values{ii});
    end
end
      
disp(outname)

% create and save the structure S
S = Z_scoord(theta_s,theta_b,tcline,hmin,N,Vtransform,Vstretching);
save(outname,'S');

% ****** Notes on S-coordinate parameters *********************
% ** THESE MUST BE IDENTICAL IN THE .IN FILE ******************
% THETA_S resolution goes to surface as thetaS increases
% THETA_B resolution to bed and surface equally at 1
% TCLINE sets the depth of the pycnocline to resove
% N number of vertical levels (rho)
% Vtransform = 1 which transformation equation to apply (1 or 2)
% Vstretching = 1 which stretching paramaterization to apply (1, 2, or 3).
% See: https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
% for transformation information.
% Also note that if hc = 0, Vtransforms 1 and 2 are identical
% *************************************************************


