function [varstruct] = Z_make_varstruct(varname)

% 10/27/2008  Parker MacCready
% This provides the structure "varstruct" required by SNCTOOLS to create a
% varible in a netcdf file.
%
% It has ROMS 3.x variable names, and tries to give them the same
% attributes as you would find in a ROMS history nc file (for grid
% varibles) or following input values from External/varinfo.dat.

varstruct.Name = varname;

switch varname
    % ********** GRID VARIABLES ***********************
    case 'spherical'
        varstruct.Nctype = 'char';
        long_name = 'grid type logical switch';
        option_T = 'spherical';
        option_F = 'Cartesian';
        varstruct.Attribute = struct('Name', ...
            {'long_name','option_T','option_F'}, ...
            'Value',{long_name,option_T,option_F});
    case 'xl'
        varstruct.Nctype = 'double';
        long_name = 'domain length in the XI-direction';
        units = 'meter';
        varstruct.Attribute = struct('Name',{'long_name','units'}, ...
            'Value',{long_name,units});
    case 'el'
        varstruct.Nctype = 'double';
        long_name = 'domain length in the ETA-direction';
        units = 'meter';
        varstruct.Attribute = struct('Name',{'long_name','units'}, ...
            'Value',{long_name,units});
    case {'h', 'hraw', 'hcarve'}
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'bathymetry at RHO-points';
        units = 'meter';
        coordinates = 'lat_rho lon_rho';
        field = 'bath, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'f'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'Coriolis parameter at RHO-points';
        units = 'second-1';
        coordinates = 'lat_rho lon_rho';
        field = 'coriolis, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'pm'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'curvilinear coordinate metric in XI';
        units = 'meter-1';
        coordinates = 'lat_rho lon_rho';
        field = 'pm, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'pn'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'curvilinear coordinate metric in ETA';
        units = 'meter-1';
        coordinates = 'lat_rho lon_rho';
        field = 'pn, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'x_rho'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'x position of RHO-points';
        units = 'meter';
        field = 'x_rho, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'y_rho'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'y position of RHO-points';
        units = 'meter';
        field = 'y_rho, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'x_u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_u';'xi_u'};
        long_name = 'x position of U-points';
        units = 'meter';
        field = 'x_u, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'y_u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_u';'xi_u'};
        long_name = 'y position of U-points';
        units = 'meter';
        field = 'y_u, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'x_v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_v';'xi_v'};
        long_name = 'x position of V-points';
        units = 'meter';
        field = 'x_v, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'y_v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_v';'xi_v'};
        long_name = 'y position of V-points';
        units = 'meter';
        field = 'y_v, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'x_psi'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_psi';'xi_psi'};
        long_name = 'x position of PSI-points';
        units = 'meter';
        field = 'x_psi, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'y_psi'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_psi';'xi_psi'};
        long_name = 'y position of PSI-points';
        units = 'meter';
        field = 'y_psi, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lon_rho'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'longitude of RHO-points';
        units = 'degree_east';
        field = 'lon_rho, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lat_rho'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'latitude of RHO-points';
        units = 'degree_north';
        field = 'lat_rho, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lon_u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_u';'xi_u'};
        long_name = 'longitude of U-points';
        units = 'degree_east';
        field = 'lon_u, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lat_u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_u';'xi_u'};
        long_name = 'latitude of U-points';
        units = 'degree_north';
        field = 'lat_u, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lon_v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_v';'xi_v'};
        long_name = 'longitude of V-points';
        units = 'degree_east';
        field = 'lon_v, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lat_v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_v';'xi_v'};
        long_name = 'latitude of V-points';
        units = 'degree_north';
        field = 'lat_v, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lon_psi'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_psi';'xi_psi'};
        long_name = 'longitude of PSI-points';
        units = 'degree_east';
        field = 'lon_psi, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'lat_psi'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_psi';'xi_psi'};
        long_name = 'latitude of PSI-points';
        units = 'degree_north';
        field = 'lat_psi, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'mask_rho'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_rho';'xi_rho'};
        long_name = 'mask on RHO-points';
        option_0 = 'land';
        option_1 = 'water';
        coordinates = 'lat_rho lon_rho';
        varstruct.Attribute = struct('Name', ...
            {'long_name','option_0','option_1','coordinates'}, ...
            'Value',{long_name,option_0,option_1,coordinates});
    case 'mask_u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_u';'xi_u'};
        long_name = 'mask on U-points';
        option_0 = 'land';
        option_1 = 'water';
        coordinates = 'lat_u lon_u';
        varstruct.Attribute = struct('Name', ...
            {'long_name','option_0','option_1','coordinates'}, ...
            'Value',{long_name,option_0,option_1,coordinates});
    case 'mask_v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_v';'xi_v'};
        long_name = 'mask on V-points';
        option_0 = 'land';
        option_1 = 'water';
        coordinates = 'lat_v lon_v';
        varstruct.Attribute = struct('Name', ...
            {'long_name','option_0','option_1','coordinates'}, ...
            'Value',{long_name,option_0,option_1,coordinates});
    case 'mask_psi'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'eta_psi';'xi_psi'};
        long_name = 'mask on psi-points';
        option_0 = 'land';
        option_1 = 'water';
        coordinates = 'lat_psi lon_psi';
        varstruct.Attribute = struct('Name', ...
            {'long_name','option_0','option_1','coordinates'}, ...
            'Value',{long_name,option_0,option_1,coordinates});
        % ********** END GRID VARIABLES ********************
        % ********** TIDAL VARIABLES ***********************
    case 'tide_period'
        varstruct.Name = 'tide_period';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period'};
        long_name = 'tidal period';
        units = 'hours';
        coordinates = 'tide_period';
        field = 'tide_period, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Eamp'
        varstruct.Name = 'tide_Eamp';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'tidal elevation amplitude';
        units = 'meter';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Eamp, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Ephase'
        varstruct.Name = 'tide_Ephase';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'tidal elevation phase angle';
        units = 'degrees';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Ephase, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Cphase'
        varstruct.Name = 'tide_Cphase';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'tidal current phase angle';
        units = 'degrees';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Cphase, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Cangle'
        varstruct.Name = 'tide_Cangle';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'tidal current inclination angle';
        units = 'degrees';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Cangle, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Cmax'
        varstruct.Name = 'tide_Cmax';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'maximum tidal current, ellipse major axis';
        units = 'meter second-1';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Cmax, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
    case 'tide_Cmin'
        varstruct.Name = 'tide_Cmin';
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'tide_period';'eta_rho';'xi_rho'};
        long_name = 'minimum tidal current, ellipse minor axis';
        units = 'meter second-1';
        coordinates = 'lat_rho lon_rho tide_period';
        field = 'tide_Cmin, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','coordinates','field'}, ...
            'Value',{long_name,units,coordinates,field});
        % ********** END TIDAL VARIABLES ********************
        % ********** CLIM VARIABLES *************************
   case 'salt'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'t','d','y','x'};
        long_name = 'Salinity';
        units = 'psu';
        field = 'S, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'temp'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'t','d','y','x'};
        long_name = 'Temperature';
        units = 'degrees-C';
        field = 'T, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'u'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'t','d','y','x'};
        long_name = 'east velocity';
        units = 'meter second-1';
        field = 'U, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'v'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'t','d','y','x'};
        long_name = 'north velocity';
        units = 'meter second-1';
        field = 'V, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field}); 
    case 'zeta'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'t','y','x'};
        long_name = 'sea surface height';
        units = 'meters';
        field = 'SSH, scalar';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});   
    case 'age_01'  %added to ocean_ini so dimensions different
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'time','s_rho','eta_rho','xi_rho'};
        long_name = 'water (age) concentration';
        units = 'nondimensional';
        field = 'age_, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});
    case 'age_02'
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'time','s_rho','eta_rho','xi_rho'};
        long_name = 'water (age) concentration';
        units = 'metr3 meter-3';
        field = 'age_, scalar, series';
        varstruct.Attribute = struct('Name', ...
            {'long_name','units','field'}, ...
            'Value',{long_name,units,field});  
        % ********** END CLIM VARIABLES **********************
end
