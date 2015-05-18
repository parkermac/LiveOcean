function clm2bry(out_dir, cname, bname)
% clm2bry.m  6/12/2007  Parker MacCready
%
% edited by DAS, 1/18/2009 to use snctools for roms 3.#
% edited by SNG, 2/16/2011 to change clmname = [out_dir, cname];
% edited by PM 11/18/2014 to not mess with time convention
%
% this creates a boundary netcdf file for ROMS using a climatology
% netcdf file
%
% based on code written by David Darr and others
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bryname = [out_dir, bname];
clmname = [out_dir, cname]; 

% get sizes
dinfo=nc_getdiminfo(clmname, 'xi_rho');    Lrho = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 'eta_rho');   Mrho = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 's_rho');     nlay = dinfo.Length;
dinfo=nc_getdiminfo(clmname, 'salt_time'); ntime = dinfo.Length;

% the naming here is a leftover from Peneven
LLm = Lrho-2;
MMm = Mrho-2;
% all the -1's below are because of 0 index starting

for ivar = 1:7
    switch ivar
        case 1
            var = 'salt'; 
            var_n = squeeze(nc_varget(clmname,var,[0 0 MMm+2-1 0],[ntime nlay 1 LLm+2]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay 1 LLm+2]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 0 LLm+2-1],[ntime nlay MMm+2 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay MMm+2 1]));
        case 2
            var = 'temp';
            var_n = squeeze(nc_varget(clmname,var,[0 0 MMm+2-1 0],[ntime nlay 1 LLm+2]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay 1 LLm+2]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 0 LLm+2-1],[ntime nlay MMm+2 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay MMm+2 1]));
        case 3
            var = 'zeta';
            var_n = squeeze(nc_varget(clmname,var,[0 MMm+2-1 0],[ntime 1 LLm+2]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0],[ntime 1 LLm+2]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 LLm+2-1],[ntime MMm+2 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0],[ntime MMm+2 1]));
        case 4
            var = 'u';
            var_n = squeeze(nc_varget(clmname,var,[0 0 MMm+2-1 0],[ntime nlay 1 LLm+1]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay 1 LLm+1]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 0 LLm+1-1],[ntime nlay MMm+2 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay MMm+2 1]));
        case 5
            var = 'v';
            var_n = squeeze(nc_varget(clmname,var,[0 0 MMm+1-1 0],[ntime nlay 1 LLm+2]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay 1 LLm+2]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 0 LLm+2-1],[ntime nlay MMm+1 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0 0],[ntime nlay MMm+1 1]));
        case 6
            var = 'ubar';
            var_n = squeeze(nc_varget(clmname,var,[0 MMm+2-1 0],[ntime 1 LLm+1]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0],[ntime 1 LLm+1]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 LLm+1-1],[ntime MMm+2 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0],[ntime MMm+2 1]));
        case 7
            var = 'vbar';
            var_n = squeeze(nc_varget(clmname,var,[0 MMm+1-1 0],[ntime 1 LLm+2]));
            var_s = squeeze(nc_varget(clmname,var,[0 0 0],[ntime 1 LLm+2]));
            var_e = squeeze(nc_varget(clmname,var,[0 0 LLm+2-1],[ntime MMm+1 1]));
            var_w = squeeze(nc_varget(clmname,var,[0 0 0],[ntime MMm+1 1]));
    end
    % rename things
    eval([var,'_north = var_n;']);
    eval([var,'_south = var_s;']);
    eval([var,'_east = var_e;']);
    eval([var,'_west = var_w;']);
end

% get some time vectors
scale = 1;
salt_time = nc_varget(clmname, 'salt_time')/scale;
temp_time = nc_varget(clmname, 'temp_time')/scale;
zeta_time = nc_varget(clmname, 'temp_time')/scale;
v3d_time = nc_varget(clmname, 'v3d_time')/scale;
v2d_time = nc_varget(clmname, 'v2d_time')/scale;

% *****************************************************************
% create the boundary file
my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
nc_create_empty( bryname, my_mode );
% global attributes
nc_padheader (bryname, 10000 );
nc_attput(bryname, nc_global, 'type','ROMS Boundary File');
nc_attput(bryname, nc_global, 'history',['from ' out_dir 'ocean_clm.nc']);

% define dimensions
nc_add_dimension(bryname, 'temp_time', length(temp_time)); 
nc_add_dimension(bryname, 'salt_time', length(salt_time));
nc_add_dimension(bryname, 'zeta_time', length(zeta_time));
nc_add_dimension(bryname, 'v3d_time', length(v3d_time));
nc_add_dimension(bryname, 'v2d_time', length(v2d_time));
%
nc_add_dimension(bryname, 's_rho', nlay);
nc_add_dimension(bryname, 'xi_rho', LLm+2);
nc_add_dimension(bryname, 'eta_rho', MMm+2);
nc_add_dimension(bryname, 'xi_u', LLm+1);
nc_add_dimension(bryname, 'eta_u', MMm+2);
nc_add_dimension(bryname, 'xi_v', LLm+2);
nc_add_dimension(bryname, 'eta_v', MMm+1);
%

% this sets the units for all the time dimensions
units = ['seconds'];
varstruct.Attribute = struct('Name', ...
    {'units'},'Value',{units});

varstruct.Name = 'salt_time'; varstruct.Dimension = {'salt_time'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'temp_time'; varstruct.Dimension = {'temp_time'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'zeta_time'; varstruct.Dimension = {'zeta_time'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v3d_time'; varstruct.Dimension = {'v3d_time'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v2d_time'; varstruct.Dimension = {'v2d_time'};
nc_addvar(bryname, varstruct);
%
varstruct.Name = 'salt_north'; varstruct.Dimension = {'salt_time','s_rho','xi_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'temp_north'; varstruct.Dimension = {'temp_time','s_rho','xi_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'u_north'; varstruct.Dimension = {'v3d_time','s_rho','xi_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v_north'; varstruct.Dimension = {'v3d_time','s_rho','xi_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'ubar_north'; varstruct.Dimension = {'v2d_time','xi_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'vbar_north'; varstruct.Dimension = {'v2d_time','xi_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'zeta_north'; varstruct.Dimension = {'zeta_time','xi_rho'};
nc_addvar(bryname, varstruct);
%
varstruct.Name = 'salt_south'; varstruct.Dimension = {'salt_time','s_rho','xi_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'temp_south'; varstruct.Dimension = {'temp_time','s_rho','xi_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'u_south'; varstruct.Dimension = {'v3d_time','s_rho','xi_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v_south'; varstruct.Dimension = {'v3d_time','s_rho','xi_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'ubar_south'; varstruct.Dimension = {'v2d_time','xi_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'vbar_south'; varstruct.Dimension = {'v2d_time','xi_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'zeta_south'; varstruct.Dimension = {'zeta_time','xi_rho'};
nc_addvar(bryname, varstruct);
%
varstruct.Name = 'salt_east'; varstruct.Dimension = {'salt_time','s_rho','eta_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'temp_east'; varstruct.Dimension = {'temp_time','s_rho','eta_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'u_east'; varstruct.Dimension = {'v3d_time','s_rho','eta_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v_east'; varstruct.Dimension = {'v3d_time','s_rho','eta_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'ubar_east'; varstruct.Dimension = {'v2d_time','eta_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'vbar_east'; varstruct.Dimension = {'v2d_time','eta_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'zeta_east'; varstruct.Dimension = {'zeta_time','eta_rho'};
nc_addvar(bryname, varstruct);
%
varstruct.Name = 'salt_west'; varstruct.Dimension = {'salt_time','s_rho','eta_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'temp_west'; varstruct.Dimension = {'temp_time','s_rho','eta_rho'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'u_west'; varstruct.Dimension = {'v3d_time','s_rho','eta_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'v_west'; varstruct.Dimension = {'v3d_time','s_rho','eta_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'ubar_west'; varstruct.Dimension = {'v2d_time','eta_u'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'vbar_west'; varstruct.Dimension = {'v2d_time','eta_v'};
nc_addvar(bryname, varstruct);
varstruct.Name = 'zeta_west'; varstruct.Dimension = {'zeta_time','eta_rho'};
nc_addvar(bryname, varstruct);

% store the data
nc_varput(bryname, 'temp_time', temp_time);
nc_varput(bryname, 'salt_time', salt_time);
nc_varput(bryname, 'zeta_time', zeta_time);
nc_varput(bryname, 'v3d_time', v3d_time);
nc_varput(bryname, 'v2d_time', v2d_time);
%
nc_varput(bryname, 'temp_north', temp_north);
nc_varput(bryname, 'salt_north', salt_north);
nc_varput(bryname, 'zeta_north', zeta_north);
nc_varput(bryname, 'v_north', v_north);
nc_varput(bryname, 'u_north', u_north);
nc_varput(bryname, 'vbar_north', vbar_north);
nc_varput(bryname, 'ubar_north', ubar_north);
%
nc_varput(bryname, 'temp_south', temp_south);
nc_varput(bryname, 'salt_south', salt_south);
nc_varput(bryname, 'zeta_south', zeta_south);
nc_varput(bryname, 'v_south', v_south);
nc_varput(bryname, 'u_south', u_south);
nc_varput(bryname, 'vbar_south', vbar_south);
nc_varput(bryname, 'ubar_south', ubar_south);
%
nc_varput(bryname, 'temp_east', temp_east);
nc_varput(bryname, 'salt_east', salt_east);
nc_varput(bryname, 'zeta_east', zeta_east);
nc_varput(bryname, 'v_east', v_east);
nc_varput(bryname, 'u_east', u_east);
nc_varput(bryname, 'vbar_east', vbar_east);
nc_varput(bryname, 'ubar_east', ubar_east);
%
nc_varput(bryname, 'temp_west', temp_west);
nc_varput(bryname, 'salt_west', salt_west);
nc_varput(bryname, 'zeta_west', zeta_west);
nc_varput(bryname, 'v_west', v_west);
nc_varput(bryname, 'u_west', u_west);
nc_varput(bryname, 'vbar_west', vbar_west);
nc_varput(bryname, 'ubar_west', ubar_west);


disp('  * DONE writing boundary file')
disp(' ')

