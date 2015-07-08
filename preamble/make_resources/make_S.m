function make_S(gridname)
% make_S.m Parker MacCready
%
% makes and saves the structure "S" for a given run

tag = 'not_needed';
addpath('../../alpha'); Ldir = Lstart(gridname, tag);

outdir = [Ldir.res,Ldir.gridname,'/'];
if ~exist(outdir,'dir'); mkdir(outdir); end;
outname = [outdir,'S.mat'];

% determine what hmin was used in the grid
gfile = [outdir,'grid.nc'];
h = nc_varget(gfile,'h');
hmin = min(h(:)); % hmin

% get the remaining S-coordinate parameters from a text file
fid = fopen([outdir,'S_COORDINATE_INFO.csv'],'r');
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


