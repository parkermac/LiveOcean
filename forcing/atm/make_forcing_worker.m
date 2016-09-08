function make_forcing_worker(gridname, tag, date_string, run_type, outdir)
%% make_forcing_worker.m
%
% ****************** for atm ***********************
%%
% I am trying a clearer workflow structure as an experiment that could work
% for any of these "reformatting" jobs that come up so often.
%
% 1. Create vectors of output times that map to vectors of input files
% 2. For each time:
%   - get the input file
%   - read in all the variables
%   - process as needed (rotate, rescale, make relative humidity)
%   - regrid to the output grid (here just horizontal)
% 3. Write each variable to its place in the netcdf output file

% Notes: need to tell this where to find the WRF files, and
% specify year, month, day

addpath('../../alpha'); Ldir = Lstart(gridname, tag);
start_time = datenum(now);

%% atm-specific code
addpath('./atm_fun');

% the forecast day
yrs = date_string(1:4);
mos = date_string(6:7);
dys = date_string(9:10);
yr = str2double(yrs);
mo = str2double(mos);
dy = str2double(dys);

td_now = datenum(yr,mo,dy);
td_prev = td_now - 1;
td_next = td_now + 1;

%% create the output time vector

if strcmp(run_type,'backfill')
    % This deals with a single case where we only have 24 hours of WRF
    % output, because we made this day by hand using overlap from
    % surrounding days. It also makes the backfill process faster in
    % general.
    hr_vec = 0:24;
elseif strcmp(run_type,'forecast')
    hr_vec = 0:72;
    % this gives 73 hourly values (three days, including endpoints)
end

%
dt_out = datenum(yr,mo,dy,hr_vec,0,0);
vartime = (dt_out - datenum(1970,1,1))*86400; % seconds since 1/1/1970

%% create the corresponding input files
%
% NOTE: for now let's just use the d2 (12 km) grid

% and get the parent
which_home = getenv('HOME');
switch which_home
    case '/Users/PM5'
        Info.wrf_dir = '/Users/PM5/Documents/LiveOcean_data/wrf/';
    case '/home/parker'
        Info.wrf_dir = '/pmr2/darr/wrf_crons/wrfout/';
    otherwise
        disp('Show me the way to get home')
end

indir00 = [Info.wrf_dir,yrs,mos,dys,'00/'];
% indir12 = [Info.wrf_dir,yrs,mos,dys,'12/']; don't use

% just use the 00 forecast
for tt = 1:length(hr_vec)
    hr = hr_vec(tt);
    hrs = ['00',num2str(hr)]; hrs = hrs(end-1:end);
    infile_list_d2{tt} = [indir00,'wrfout.ocean_d2.',yrs,mos,dys, ...
        '00.f',hrs,'.0000'];
    forecast_hour{tt} = hrs; % used for rain calculation
end

%% variable name lists

invar_list = {'PSFC','RAINC','RAINNC','SWDOWN','GLW', ...
    'T2','Q2','U10','V10'};

%% grids

% the ROMS grid
%gdir = [Ldir.res,Ldir.gridname,'/'];
gdir = [Ldir.data,'grids/',Ldir.gridname,'/'];
fng = [gdir,'grid.nc'];
lon = nc_varget(fng,'lon_rho');
lat = nc_varget(fng,'lat_rho');

% the WRF grids (2 for d2, and 3 for d3)
if exist(infile_list_d2{1},'file')
    lon2 = double(nc_varget(infile_list_d2{1},'XLONG'));
    lat2 = double(nc_varget(infile_list_d2{1},'XLAT'));
else
    disp('Missing WRF input files!');
end

% make blank results arrays
[NR,NC] = size(lon);
[NR2,NC2] = size(lon2);
NT = length(hr_vec);
nmat = NaN * ones(NT,NR,NC); % sized for output
nmat2 = NaN * ones(NT,NR2,NC2); % sized for input

%% read input files into separate variable files

% start with all variables at one time, and end up with
% one variable at all times

for tt = 1:NT
    fn2 = infile_list_d2{tt};
    for vv = 1:length(invar_list)
        VR = invar_list{vv};
        if tt == 1; eval([VR,'2 = nmat2;']); end;
        eval([VR,'2(tt,:,:) = nc_varget(fn2,VR);']);
    end
end

%% processing

% convert Pa to mbar
PSFC2 = PSFC2/100;

% rain
%
% WRF reports the accumulation since the start of the
% forecast, so we can estimate the precipitation RATE by
% dividing by the time since the start of the forecast.
% The 0.1 gets us from mm/sec to cm/sec
RAIN2 = nmat2; % precipitation rate [cm s-1]
for tt = 1:NT
    dts =3600 * str2num(forecast_hour{tt});
    if dts == 0
        dts = 3600;
        RAIN2(tt,:,:) = 0.1 * (RAINC2(tt+1,:,:) + RAINNC2(tt+1,:,:))/dts;
    else
        RAIN2(tt,:,:) = 0.1 * (RAINC2(tt,:,:) + RAINNC2(tt,:,:))/dts;
    end
end
% convert [cm s-1] to [kg m-2 s-1]
RAIN2 = RAIN2 * 10;

% convert K to C
T22 = T22 - 273.15;

% account for reflection
SWDOWN2 = (1 - 0.1446)*SWDOWN2;

% calculate relative humidity
RH2 = Z_wmo_RH(PSFC2,T22,Q22); % relative humidity [%]

% rotate velocity to E-W and N-S
theta = NaN * lon2;
for rr = 1:NR2
    [dist,theta(rr,1:end-1)] = sw_dist(lat2(rr,:),lon2(rr,:),'km');
end
theta(:,end) = theta(:,end-1); % pad the end
alpha = - pi*theta/180; % alpha is the angle (radians) to rotate the
% coordinate system to go from the WRF grid to East-North
sa = sin(alpha); ca = cos(alpha);
ssa = nmat2; cca = nmat2;
for tt = 1:NT
    ssa(tt,:,:) = sa;
    cca(tt,:,:) = ca;
end
UR102 = U102.*cca + V102.*ssa;
VR102 = V102.*cca - U102.*ssa;

%% regrid to the model grid

invar_list = {'PSFC','RAIN','SWDOWN','GLW', ...
    'T2','RH','UR10','VR10'};

outvar_list = {'Pair','rain','swrad','lwrad_down', ...
    'Tair','Qair','Uwind','Vwind'};

for tt = 1:NT
    for vv = 1:length(invar_list)
        VR = invar_list{vv};
        vr = outvar_list{vv};
        if tt == 1; eval([vr,' = nmat;']); end;
        eval(['OO = squeeze(',VR,'2(tt,:,:));']);
        %F = TriScatteredInterp(lon2(:),lat2(:),OO(:));
        F = scatteredInterpolant(lon2(:),lat2(:),OO(:));
        eval([vr,'(tt,:,:) = F(lon,lat);']);
    end
end

%% write out to NetCDF for ROMS

for vv = 1:length(outvar_list);
    
    ncvarname = outvar_list{vv};
    
    % get grid size
    [M,L] = size(lon);
    
    % get some attributes
    [nclongname,ncunits,nctimename] = atm_attributes(ncvarname);
    
    % create the NetCDF file
    frcname = [outdir, ncvarname, '.nc'];
    
    my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
    nc_create_empty( frcname, my_mode );
    disp(' ')
    disp(['  - creating ',frcname]);
    % global attributes
    nc_padheader (frcname, 20000 );
    nc_attput(frcname, nc_global, 'type','ROMS Surface Forcing File');
    
    % define dimensions
    nc_add_dimension(frcname, 'xi_rho', L);
    nc_add_dimension(frcname, 'eta_rho', M);
    nc_add_dimension(frcname, nctimename, length(vartime));
    
    % define variables and attributes
    %   first time
    varstruct.Name = nctimename;
    varstruct.Dimension = {nctimename};
    long_name = [nclongname,' time'];
    units = 'seconds';
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    %   then the field
    varstruct.Name = ncvarname;
    varstruct.Dimension = {nctimename, 'eta_rho', 'xi_rho'};
    long_name = nclongname;
    units = ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    % write the time
    nc_varput(frcname, nctimename, vartime);
    
    % write the variable
    eval(['nc_varput(frcname, ncvarname, ',ncvarname,');']);
    disp(['  -- done filling ',frcname]);
    
end % end of variable loop

%% things to use for checking result
t_datenum = dt_out;

%% Final output
datestr_format = 'yyyy.mm.dd HH:MM:SS';
end_time = datenum(now);
fid = fopen([outdir,'Info/process_status.csv'],'w');
fprintf(fid,'%s\n',['start_time,',datestr(start_time, datestr_format)]);
fprintf(fid,'%s\n',['end_time,',datestr(end_time, datestr_format)]);
% test for existence of output files
all_files = true;
for vv = 1:length(outvar_list)
    ncvarname = outvar_list{vv};
    frcname = [outdir, ncvarname, '.nc'];
    if ~exist(frcname,'file')
        all_files = false;
    end
end
fprintf(fid,'%s\n',['var_start_time,',datestr(t_datenum(1), datestr_format)]);
fprintf(fid,'%s\n',['var_end_time,',datestr(t_datenum(end), datestr_format)]);
if all_files
    fprintf(fid,'%s\n','result,success');
else
    fprintf(fid,'%s\n','result,fail');
    fprintf(fid,'%s\n','reason,not all output files made');
end
fclose(fid);

