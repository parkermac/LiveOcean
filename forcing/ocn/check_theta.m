% check_theta.m  9/16/2015  Parker MacCready
%
% code to check on the conversion of in-situ temperature from
% HYCOM to potential temperature to use in ROMS.
%
% RESULT: it appears to work correctly, although as expected the difference
% is small, e.g.
% * theta is about 0.15  degC lower at 2 degC
% * theta is about 0.001 degC lower at 10 degC

clear
gridname = 'cascadia1';
tag = 'base';
addpath('../../alpha'); Ldir = Lstart(gridname, tag);
date_string = '2013.03.30';

indir0 = [Ldir.out,Ldir.gtag,'/f',date_string,'_ORIG/ocn/'];
indir1 = [Ldir.out,Ldir.gtag,'/f',date_string,'/ocn/'];

f0 = [indir0,'ocean_clm.nc'];
f1 = [indir1,'ocean_clm.nc'];

% REMEMBER this is zero-based counting
t0 = squeeze(nc_varget(f0,'temp',[0 0 0 0],[1 -1 -1 -1]));
t1 = squeeze(nc_varget(f1,'temp',[0 0 0 0],[1 -1 -1 -1]));

close all
figure
plot(t0(:),t1(:),'.r')
hold on
plot([0 12],[0 12],'-k')
xlabel('Original Temperature')
ylabel('Potential Temperature')
grid on
axis square

