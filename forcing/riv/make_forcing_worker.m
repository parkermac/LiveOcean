function make_forcing_worker(gridname, tag, date_string, run_type, outdir)
%% make_forcing_worker.m
%
% ****************** for riv ***********************
%
% makes rivers.nc using data created by make_forcing_main.py

addpath('../../alpha'); Ldir = Lstart(gridname, tag);
start_time = datenum(now);

%% river-specific code
addpath('./riv_fun');
% define and load preamble files
indir = [Ldir.res,Ldir.gridname,'/'];
gridfile = [indir,'grid.nc'];
load([indir,'S.mat']); % get structure "S"

% load river location Info
river_index_file = [indir,'river_indices.mat'];
load(river_index_file); % load structure "rout"
% how many rivers
for ii = 1:length(rout); uid(ii) = rout(ii).id; end;

% initialize is a flag to indicate whether you want the river flow to ramp
% up (1) or not (0) (only needed when starting a run from climatology)
initialize = 0;

% name output file
riv_file_out=[outdir,'rivers.nc'];

% input data files from main
flowFile = [outdir,'Data/current_river_table.csv'];
qraw = importdata(flowFile);
tempFile = [outdir,'Data/current_triver_table.csv'];
traw = importdata(tempFile);
% We appear to be dropping the first column (time)
Qr_flow = qraw.data(:,2:end);
T_riv = traw.data(:,2:end);
rname = qraw.colheaders(2:end); % a cell array

% we will always tell model time as seconds since 1/1/1970
river_time = qraw.data(:,1); 

% get some grid Info
h=nc_varget(gridfile,'h');
[MP,LP]=size(h);

% pass data to function that makes river NetCDF forcing file
make_river_netcdf(MP, LP, S, rout, riv_file_out, river_time, ...
    Qr_flow, T_riv, uid, initialize);

% add more tracers

%% things to use for checking result 
outvar_list = {'rivers'};
t_datenum = river_time/86400 + datenum(1970,1,1);

%% Final output
end_time = datenum(now);
fid = fopen([outdir,'Info/process_status.csv'],'w');
fprintf(fid,'%s\n',['start_time,',datestr(start_time)]);
fprintf(fid,'%s\n',['end_time,',datestr(end_time)]);
% test for existence of output files
all_files = true;
for vv = 1:length(outvar_list)
    ncvarname = outvar_list{vv};
    frcname = [outdir, ncvarname, '.nc'];
    if ~exist(frcname,'file')
        all_files = false;
    end
end
fprintf(fid,'%s\n',['var_start_time,',datestr(t_datenum(1))]);
fprintf(fid,'%s\n',['var_end_time,',datestr(t_datenum(end))]);
if all_files
    fprintf(fid,'%s\n','result,success');
else
    fprintf(fid,'%s\n','result,fail');
    fprintf(fid,'%s\n','reason,not all output files made');
end
fclose(fid);





