function make_forcing_worker(gridname, tag, date_string, run_type, outdir)
%% make_forcing_worker.m
%
% ****************** for bio ***********************
%
% makes [].nc using data created by make_forcing_main.py

%% Debugging
% gridname = 'cascadia1';
% tag = 'base';
% date_string = '2014.03.11';
% run_type = 'backfill';
% outdir = '~/Desktop/';

%%
addpath('../../alpha'); Ldir = Lstart(gridname, tag);
start_time = datenum(now);

%% bio-specific code

addpath('./bio_fun');

% define and load preamble files
%indir = [Ldir.res,Ldir.gridname,'/'];
gdir = [Ldir.data,'grids/',Ldir.gridname,'/'];
grdname = [gdir,'grid.nc'];

% define locations of existing ocn and riv files
ocn_dir = 'ocn1';
clmname = [Ldir.LOo,Ldir.gtag,'/f',date_string,'/',ocn_dir,'/ocean_clm_bio.nc'];
bryname = [Ldir.LOo,Ldir.gtag,'/f',date_string,'/',ocn_dir,'/ocean_bry_bio.nc'];
ininame = [Ldir.LOo,Ldir.gtag,'/f',date_string,'/',ocn_dir,'/ocean_ini_bio.nc'];
rivname = [Ldir.LOo,Ldir.gtag,'/f',date_string,'/riv/rivers_bio.nc'];

%%
NO3_method='PL_Salt';

%make_bio_BCsfast_withDIC(out_dir, cname, bname, rname, NO3_method, iname);
make_bio_BCsfast_withDIC(grdname, clmname, bryname, ininame, rivname, NO3_method);


%% things to use for checking result 
outfile_list = {clmname,bryname,ininame,rivname};
t_datenum = nc_varget(rivname,'river_time')/86400 + datenum(1970,1,1);

%% Final output
datestr_format = 'yyyy.mm.dd HH:MM:SS';
end_time = datenum(now);
fid = fopen([outdir,'Info/process_status.csv'],'w');
fprintf(fid,'%s\n',['start_time,',datestr(start_time, datestr_format)]);
fprintf(fid,'%s\n',['end_time,',datestr(end_time, datestr_format)]);
% test for existence of output files
all_files = true;
for vv = 1:length(outfile_list)
    frcname = outfile_list{vv};
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





