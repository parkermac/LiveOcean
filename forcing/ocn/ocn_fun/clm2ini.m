function clm2ini(out_dir,clmname1,iname,date_string)
% clm2ini.m  12/6/2011  Parker MacCready, DAS, KAD, SNG
%
% This creates a NetCDF file to be used for ROMS initialization, working
% with an existing climatology file
%
% edited 11/18/2014 by PM to net mess with time units
%------------------------------------------------------------------------

ininame = [out_dir, iname];
clmname = [out_dir, clmname1];

my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
nc_create_empty( ininame, my_mode );
% global attributes
nc_padheader (ininame, 20000 );
nc_attput(ininame, nc_global, 'type','ROMS Initial Condition File');
nc_attput(ininame, nc_global, 'history',['from ' out_dir clmname]);

%% the forecast day
yrs = date_string(1:4);
mos = date_string(6:7);
dys = date_string(9:10);
yr = str2double(yrs);
mo = str2double(mos);
dy = str2double(dys);

tini = (datenum(yr,mo,dy) - datenum(1970,1,1))* 86400;

%% assume that all climatology variable have the same time vector
% AND that the run initialization time exists in the climatology file
% tin = nc_varget(clmname,'temp_time');
itin = 1;
% tini = tin(itin);

% get sizes
dinfo=nc_getdiminfo(clmname, 'xi_rho');    Lrho = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 'eta_rho');   Mrho = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 's_rho');     nlay = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 'salt_time'); ntime = dinfo.Length;
% the naming here is a leftover from Peneven
LLm = Lrho-2;
MMm = Mrho-2;
% add dimensions to ini file
nc_add_dimension(ininame, 's_rho', nlay);
nc_add_dimension(ininame, 'time', 1);
nc_add_dimension(ininame, 'xi_rho', LLm+2);
nc_add_dimension(ininame, 'eta_rho', MMm+2);
nc_add_dimension(ininame, 'xi_u', LLm+1);
nc_add_dimension(ininame, 'eta_u', MMm+2);
nc_add_dimension(ininame, 'xi_v', LLm+2);
nc_add_dimension(ininame, 'eta_v', MMm+1);
nc_add_dimension(ininame, 'one', 1);

for ii = 1:7
    switch ii
        case 1
            varname = 'salt';
            varval(1,:,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0 0], [1 nlay Mrho Lrho]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','s_rho','eta_rho','xi_rho'};
            nc_addvar(ininame, varstruct);
        case 2
            varname = 'temp';
            varval(1,:,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0 0], [1 nlay Mrho Lrho]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','s_rho','eta_rho','xi_rho'};
            nc_addvar(ininame, varstruct);
        case 3
            varname = 'u';
            varval(1,:,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0 0], [1 nlay MMm+2 LLm+1]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','s_rho','eta_u','xi_u'};
            nc_addvar(ininame, varstruct);
        case 4
            varname = 'v';
            varval(1,:,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0 0], [1 nlay MMm+1 LLm+2]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','s_rho','eta_v','xi_v'};
            nc_addvar(ininame, varstruct);
        case 5
            varname = 'ubar';
            varval(1,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0], [1 MMm+2 LLm+1]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','eta_u','xi_u'};
            nc_addvar(ininame, varstruct);
        case 6
            varname = 'vbar';
            varval(1,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0], [1 MMm+1 LLm+2]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','eta_v','xi_v'};
            nc_addvar(ininame, varstruct);
        case 7
            varname = 'zeta';
            varval(1,:,:) = squeeze(nc_varget(clmname, ... %DAS added ()
                varname, [itin-1 0 0], [1 Mrho Lrho]));
            varstruct.Name = varname;
            varstruct.Dimension = {'time','eta_rho','xi_rho'};
            nc_addvar(ininame, varstruct);
    end
    nc_varput(ininame, varname, varval);
    clear varval
end

% create variables
varstruct.Name = 'ocean_time';
varstruct.Dimension = {'time'};
long_name = ['Initial condition time'];
units = ['seconds'];
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(ininame, varstruct);
nc_varput(ininame, 'ocean_time', tini);

varstruct.Name = 'tstart';
varstruct.Dimension = {'one'};
long_name = ['start time'];
units = ['seconds'];
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(ininame, varstruct);
nc_varput(ininame, 'tstart', tini);

varstruct.Name = 'tend';
varstruct.Dimension = {'one'};
long_name = ['end time'];
units = ['seconds'];
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(ininame, varstruct);
nc_varput(ininame, 'tend', tini);

disp('  * DONE writing initial file')
disp(' ')

