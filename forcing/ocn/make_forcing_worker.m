function make_forcing_worker(gridname, tag, date_string, run_type, outdir)
%% make_forcing_worker.m
%
% ****************** for ocn ***********************
%
% This creates a NetCDF file to be used for ROMS forcing.  It is designed
% to read in the files created by make_forcing_main.py
%
% It makes: a ROMS input climatology file
% with Information on T, S, u, v, ubar, vbar and zeta

addpath('../../alpha'); Ldir = Lstart(gridname, tag);
start_time = datenum(now);

%% ocn-specific code
addpath('./ocn_fun');

% the ROMS grid Info
%resdir = [Ldir.res,Ldir.gridname,'/'];
gdir = [Ldir.data,'grids/',Ldir.gridname,'/'];
gridfile = [gdir,'grid.nc'];
load([gdir,'S.mat']);

% this is where the processed HYCOM files are - from python
indir0 = [outdir,'Data/'];

% output name
clmname = 'ocean_clm.nc'

%% starting the work

outfile = [outdir,clmname];
clim_netcdf_new(gridfile,outfile,S.N);

h = nc_varget(gridfile,'h');
lon_rho = nc_varget(gridfile,'lon_rho');
lat_rho = nc_varget(gridfile,'lat_rho');

varname_list={'s3d';'t3d';'u3d';'v3d'}; % all except ssh

for ivv = 1:length(varname_list) % start of VARNAME LOOP
    tstart = tic;
    varname = varname_list{ivv};
    
    ts_to_get_uv = nc_varget([indir0,varname,'.nc'],'dt');
    nts_uv = length(ts_to_get_uv);
    
    [V,AAlon,AAlat,AAmask] = ...
        clim_netcdf(indir0,varname,gridfile);
    [M,L] = size(AAlon);
    hh = interp2(lon_rho,lat_rho,h,AAlon,AAlat);
    disp(['    - writing ',V.invarname,' to ',V.ncvarname]);
    
    clim_netcdf_time(outfile,V,nts_uv);
    if V.do_addtimevar
        nc_varput(outfile, V.nctimename, ts_to_get_uv);
    end
    
    switch varname
        case 'u3d'
            ubar_mat = NaN * ones(nts_uv,M,L);
        case 'v3d';
            vbar_mat = NaN * ones(nts_uv,M,L);
    end
    
    varstruct.Name = V.ncvarname;
    varstruct.Dimension = {V.nctimename, 's_rho', ...
        V.ncetaname, V.ncxiname};
    long_name = V.nclongname;
    units = V.ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(outfile, varstruct);
    fn = [indir0,varname,'.nc'];
    inz = nc_varget(fn, 'z');
    AN = length(inz);
    
    for tt = 1:nts_uv % start of TIME LOOP
        fn = [indir0,varname,'.nc'];
        if tt == 1
            % HYCOM grid
            Alon = nc_varget(fn, 'lon');
            Alon = Alon(1,:);
            Alat = nc_varget(fn, 'lat');
            Alat = Alat(:,1);
            % 9/16/2015 make a 3d z grid for HYCOM
            % to compute pressure for potential temp
            disp([' -- saving Apress for time ',num2str(tt)])
            Az3 = repmat(-inz(:),1,length(Alat),length(Alon));
            Alat_mean = nanmean(Alat(:));
            Apress = sw_pres(Az3,Alat_mean);
            % end of Az3 code
            [ALON,ALAT] = meshgrid(Alon,Alat);
            % ROMS grid
            [z_rho,z_w] = Z_s2z(hh,0*hh,S);
            DZ = diff(z_w,1,1);
        end
        
        % reinitialize the storage arrays
        % AA is the empty matrix on the ROMS grid in which to put
        % the field from this time step, BUT it has the vertical grid
        % of the HYCOM file
        AA = NaN * ones(AN,M,L);
        % AAA is the empty matrix on the ROMS grid in which to put
        % the field from this time step, and it has the vertical grid
        % of the ROMS model.
        AAA = NaN * ones(S.N,M,L);
        
        % load the variable field at this time step
        clear A
        % INTERPOLATE HORIZONTALLY TO ROMS GRIDS (fills whole domain,
        % if there is ANY data on that level)
        A = nc_varget(fn, V.invarname,[tt-1 0 0 0],...
            [1 -1 -1 -1]);
        % 9/16/2015 Change in situ temperature to potential temperature
        switch varname
            case 's3d'
                disp([' -- saving Asalt for time ',num2str(tt)])
                Asalt(tt,:,:,:) = A; % to use in sw_ptmp call below
                % we are assuming that s3d comes before t3d, as it should
                % (should preallocate)
            case 't3d'
                disp([' -- calculating ptmp for time ',num2str(tt)])
                Aorig = A;
                A = sw_ptmp(squeeze(Asalt(tt,:,:,:)),Aorig,Apress,0*Apress);
        end
        
        for zz = 1:AN % start of Z-LOOP
            this_A = squeeze(A(zz,:,:));
            AA(zz,:,:) = interp2(ALON,ALAT,this_A,AAlon,AAlat);
        end
        
        % INTERPOLATE VERTICALLY TO ROMS S-SURFACES
        fillv = 0;
        this_z_ncom = inz;
        for jj = 1:M % start of XY-LOOP
            for ii = 1:L
                % only do interpolation on non-masked points
                if AAmask(jj,ii)
                    this_AA = AA(:,jj,ii);
                    this_z_roms = z_rho(:,jj,ii);
                    this_AAA = interp1q(this_z_ncom, this_AA, ...
                        this_z_roms);
                    AAA(:,jj,ii) = this_AAA;
                else
                    AAA(:,jj,ii) = fillv;
                end
            end
        end
        
        switch varname
            case 'u3d'
                ubar_mat(tt,:,:) = squeeze(sum(AAA.*DZ,1)) ./ hh;
                if tt == nts_uv
                    tempufile = [outdir 'temporary_ubar_storage.mat'];
                    save(tempufile,'ubar_mat')
                end
            case 'v3d'
                vbar_mat(tt,:,:) = squeeze(sum(AAA.*DZ,1)) ./ hh;
                if tt == nts_uv
                    tempvfile = [outdir 'temporary_vbar_storage.mat'];
                    save(tempvfile,'vbar_mat')
                end
        end
        
        % write the field at one time level
        [nz,nx,ny]=size(AAA);
        nc_varput(outfile, V.ncvarname, AAA, [tt-1 0 0 0], [1 nz nx ny]);
        clear nz nx ny AAA
        
    end % end of TIME LOOP
    delta_t = toc(tstart);
    disp(['      ',num2str(round(delta_t)),' sec for ', ...
        num2str(nts_uv),' times'])
end

% ************************ SSH **********************************
varname = 'ssh';

ts_to_get = nc_varget([indir0,varname,'.nc'],'dt');
nts = length(ts_to_get);

[V,AAlon,AAlat,AAmask] = ...
    clim_netcdf(indir0,varname,gridfile);
[M,L] = size(AAlon);
disp(['    - writing ',V.invarname,' to ',V.ncvarname]);

clim_netcdf_time(outfile,V,nts);
nc_varput(outfile, V.nctimename, ts_to_get);

varstruct.Name = V.ncvarname;
varstruct.Dimension = {V.nctimename, V.ncetaname, V.ncxiname};
long_name = V.nclongname;
units = V.ncunits;
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(outfile, varstruct);

for tt = 1:nts % start of TIME LOOP
    fn = [indir0,varname,'.nc'];
    
    if tt == 1
        % these are the NCOM lat and lon for this variable
        Alon = nc_varget(fn, 'lon');
        Alon = Alon(1,:);
        Alat = nc_varget(fn, 'lat');
        Alat = Alat(:,1);
        [ALON,ALAT] = meshgrid(Alon,Alat);
    end
    
    % reinitialize the storage arrays
    AA = NaN * ones(M,L); % ROMS grid
    % load the variable field at this time step
    clear A
    % INTERPOLATE HORIZONTALLY TO ROMS GRIDS (fills whole domain)
    A = nc_varget(fn, V.invarname,[tt-1 0 0],[1 -1 -1]);
    AA = interp2(ALON,ALAT,A,AAlon,AAlat);
    %fill any missing values remaining with 0
    AA(isnan(AA)) = 0;
    % write the field at one time level
    [nx,ny]=size(AA);
    nc_varput(outfile, V.ncvarname, AA, [tt-1 0 0], [1 nx ny]);
    clear nx ny AA
    
end % end of TIME LOOP

disp(['      ',num2str(nts),' times'])
% ************************ END SSH ******************************

varname_list={'ubar';'vbar'};
%varname_list={'ubar'};%;'vbar'};
for ivv = 1:length(varname_list) % start of VARNAME LOOP
    varname = varname_list{ivv};
    [V,AAlon,AAlat,AAmask] = ...
        clim_netcdf(indir0,varname,gridfile);
    [M,L] = size(AAlon);
    disp(['    - writing ',V.invarname,' to ',V.ncvarname]);
    
    clim_netcdf_time(outfile,V,nts_uv);
    if V.do_addtimevar
        nc_varput(outfile, V.nctimename, ts_to_get_uv);
    end
    
    varstruct.Name = V.ncvarname;
    varstruct.Dimension = {V.nctimename, V.ncetaname, V.ncxiname};
    long_name = V.nclongname;
    units = V.ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(outfile, varstruct);
    
    switch varname
        case 'ubar'
            load(tempufile)
            nc_varput(outfile,V.ncvarname,ubar_mat);
            delete(tempufile);
        case 'vbar'
            load(tempvfile)
            nc_varput(outfile,V.ncvarname,vbar_mat);
            delete(tempvfile);
    end
    
end

disp('  * DONE writing climatology file')
disp(' ')

%% make other files

clm2ini(outdir,clmname,'ocean_ini.nc',date_string);
clm2bry(outdir,clmname,'ocean_bry.nc');

% add more tracers

%% things to use for checking result
outvar_list = {'ocean_clm','ocean_ini','ocean_bry'};
t_datenum = ts_to_get/86400 + datenum(1970,1,1);

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

