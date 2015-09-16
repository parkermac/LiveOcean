function make_bio_BCsfast_withDIC(out_dir, cname, bname, rname, NO3_method, iname)

%
% make_bio_BCs.m  12/22/2011 Samantha Siedlecki
%
% This code edits the ROMS forcing files (clm, bry, ini, and river) created
% previously by make_clim.m, clm2bry.m, and clm2ini.m, and make_rivers.m
% to include biological variables.
%
% Current Issues:  Please READ!!
% If the clm/bry/ini files already have biovars, then you get an error when
% make_bio_BCs tried to create new bio variables with the same name as
% existing bio vars, so you must start with clm/bry/ini files that do not
% have the bio vars added yet.

% Updated 4/25/2011 by KAD to write bio variables to rivers.nc
% out_dir: directory for output to go - for example:
%  out_dir='/home/disk/edward/siedlesa/turbo40forcingfiles/2007/';
% cname : name of the clim file - for example, cname='ocean_clm_1_test.nc';
% bname: name of the bry file - for example, bname='ocean_bry_1_test.nc';
% rname: name of the river file- for example, rname='rivers.nc';
% NO3_method: always NO3_method='PL_Salt';
% iname: name of ini file - for example, iname='ocean_ini_test.nc';
% gname: name of the grid file  - needs to be in the out_dir too.

% Updated 5/27/2015 by PM to fit into LiveOcean framework

%-----------------------------------------------------------------------

% name of climatology file
clmname = [out_dir, cname];
gname=[out_dir,'grid.nc'];
%
% check to see if this file exists already, if not display warning and exit
if ~exist(clmname,'file')
    disp(['WARNING: Climatology file does not exist.']);
end

% If the clm file does exist then proceed
if exist(clmname,'file')
    tic
    
    %open the netcdf file (assigns an id to it to access it)
    fileattrib(clmname,'+w');
    ncid=netcdf.open(clmname, 'NC_WRITE');
    %set the fill value to nothing so it will be zeros by default
    netcdf.setFill(ncid,'NC_NOFILL');
    
    %get id for variables of interest
    salt_time_id = netcdf.inqVarID(ncid,'salt_time');
    salt_id = netcdf.inqVarID(ncid,'salt');
    temp_id = netcdf.inqVarID(ncid,'temp');
    %get the info about the salt including dimensions to copy to tracers
    [~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_id);
    
    %Create time fields
    vartime = netcdf.getVar(ncid,salt_time_id,'double'); % basing bio time field off salt_time
    clm_salt = netcdf.getVar(ncid,salt_id,'double');
    clm_temp = netcdf.getVar(ncid,temp_id,'double');
    
    % Create NO3 field
    if strcmp(NO3_method,'PL_Salt')
        % Calculate NO3
        NO3 = make_NO3_field(NO3_method,clm_salt);
        % Calculate oxygen
        oxygen = make_oxy_field(NO3_method,clm_salt);
        [DIC, Talk] = make_DIC_field(NO3_method,clm_salt,clm_temp);
        
    end
    NO3=-1000*ones(size(clm_salt)); % what? PM
    
    % Create phytoplankton field
    phytoplankton = 0.01*ones(size(clm_salt)); % "Seed stock" of 0.01 microM of N
    
    % Create zooplankton field
    zooplankton = 0.01*ones(size(clm_salt)); % "Seed stock" of 0.01 microM of N
    
    % Create detritus field
    detritus = zeros(size(clm_salt)); % Initializing with no detritus
    %
    % Create detritus field
    Ldetritus = zeros(size(clm_salt)); % Initializing with no detritus
    PIC = zeros(size(clm_salt)); % Initializing with no detritus
    
    %%% WRITE VARIABLES TO NETCDF FILE %%%%%%%%%%%%%%%%%%%%%%
    %%Define new time dimensions in clm file
    
    netcdf.reDef(ncid);
    dimidn=netcdf.defDim(ncid, 'NO3_time', length(vartime));
    dimidp=netcdf.defDim(ncid, 'phyt_time', length(vartime));
    dimidz=netcdf.defDim(ncid, 'zoop_time', length(vartime));
    dimidd=netcdf.defDim(ncid, 'detritus_time', length(vartime));
    dimidLd=netcdf.defDim(ncid, 'Ldetritus_time', length(vartime));
    dimidPIC=netcdf.defDim(ncid, 'CaCO3_time', length(vartime));
    dimido=netcdf.defDim(ncid, 'oxygen_time', length(vartime));
    dimida=netcdf.defDim(ncid, 'talk_time', length(vartime));
    dimiddic=netcdf.defDim(ncid, 'tic_time', length(vartime));
    
    
    %Define the variables based on salt dimensions
    varidnt=netcdf.defVar(ncid, 'NO3_time', 'double', dimidn);
    varidpt=netcdf.defVar(ncid, 'phyt_time', 'double', dimidp);
    varidzt=netcdf.defVar(ncid, 'zoop_time', 'double', dimidz);
    variddt=netcdf.defVar(ncid, 'detritus_time', 'double', dimidd);
    varidLdt=netcdf.defVar(ncid, 'Ldetritus_time', 'double', dimidLd);
    varidPICt=netcdf.defVar(ncid, 'CaCO3_time', 'double', dimidPIC);
    varidot=netcdf.defVar(ncid, 'oxygen_time', 'double', dimido);
    varidat=netcdf.defVar(ncid, 'talk_time', 'double', dimida);
    
    variddict=netcdf.defVar(ncid, 'dic_time', 'double', dimiddic);
    
    varidn=netcdf.defVar(ncid, 'NO3',xtype,dimids);
    varidp=netcdf.defVar(ncid, 'phytoplankton', xtype,dimids);
    varidz=netcdf.defVar(ncid, 'zooplankton',  xtype,dimids);
    varidd=netcdf.defVar(ncid, 'detritus',  xtype,dimids);
    varidLd=netcdf.defVar(ncid, 'Ldetritus',  xtype,dimids);
    varidLd=netcdf.defVar(ncid, 'Ldetritus',  xtype,dimids);
    varidPIC=netcdf.defVar(ncid, 'CaCO3',  xtype,dimids);
    varido=netcdf.defVar(ncid, 'oxygen',  xtype,dimids);
    varida=netcdf.defVar(ncid, 'alkalinity',  xtype,dimids);
    variddic=netcdf.defVar(ncid, 'TIC',  xtype,dimids);
    
    
    %add the desired attributes
    netcdf.putAtt(ncid,varidnt,'long_name','NO3 time');
    netcdf.putAtt(ncid,varidnt,'units','s');
    netcdf.putAtt(ncid,varidpt,'long_name','Phytoplankton time');
    netcdf.putAtt(ncid,varidpt,'units','s');
    netcdf.putAtt(ncid,varidzt,'long_name','Zooplankton time');
    netcdf.putAtt(ncid,varidzt,'units','s');
    netcdf.putAtt(ncid,variddt,'long_name','Detritus time');
    netcdf.putAtt(ncid,variddt,'units','s');
    netcdf.putAtt(ncid,varidLdt,'long_name','Large Detritus time');
    netcdf.putAtt(ncid,varidLdt,'units','s');
    netcdf.putAtt(ncid,varidPICt,'long_name','CaCO3 time');
    netcdf.putAtt(ncid,varidPICt,'units','s');
    netcdf.putAtt(ncid,varidot,'long_name','Oxygen time');
    netcdf.putAtt(ncid,varidot,'units','s');
    netcdf.putAtt(ncid,varidat,'long_name','alkalinity time');
    netcdf.putAtt(ncid,varidat,'units','s');
    netcdf.putAtt(ncid,variddict,'long_name','TIC time');
    netcdf.putAtt(ncid,variddict,'units','s');
    
    
    netcdf.putAtt(ncid,varidn,'long_name','NO3 climatology');
    netcdf.putAtt(ncid,varidn,'units','MicroMolar');
    netcdf.putAtt(ncid,varidp,'long_name','Phytoplankton climatology');
    netcdf.putAtt(ncid,varidp,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidz,'long_name','Zooplankton climatology');
    netcdf.putAtt(ncid,varidz,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidd,'long_name','Detritus climatology');
    netcdf.putAtt(ncid,varidd,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidLd,'long_name','Large Detritus climatology');
    netcdf.putAtt(ncid,varidLd,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidPIC,'long_name','CaCO3 climatology');
    netcdf.putAtt(ncid,varidPIC,'units','MicroMolar C');
    netcdf.putAtt(ncid,varido,'long_name','Oxygen climatology');
    netcdf.putAtt(ncid,varido,'units','MicroMolar O');
    netcdf.putAtt(ncid,varida,'long_name','Talk climatology');
    netcdf.putAtt(ncid,varida,'units','MicroMolar ');
    netcdf.putAtt(ncid,variddic,'long_name','DIC climatology');
    netcdf.putAtt(ncid,variddic,'units','MicroMolar C');
    
    %%finish defining
    netcdf.endDef(ncid);
    
    %% Store the data
    netcdf.putVar(ncid,varidnt, vartime);
    netcdf.putVar(ncid,varidpt, vartime);
    netcdf.putVar(ncid,varidzt, vartime);
    netcdf.putVar(ncid,variddt, vartime);
    netcdf.putVar(ncid,varidLdt, vartime);
    netcdf.putVar(ncid,varidPICt, vartime);
    netcdf.putVar(ncid,varidot, vartime);
    netcdf.putVar(ncid,varidat, vartime);
    netcdf.putVar(ncid,variddict, vartime);
    
    netcdf.putVar(ncid, varidn, NO3);
    netcdf.putVar(ncid, varidp, phytoplankton);
    netcdf.putVar(ncid, varidz, zooplankton);
    netcdf.putVar(ncid, varidd, detritus);
    netcdf.putVar(ncid, varidLd, Ldetritus);
    netcdf.putVar(ncid, varidPIC, PIC);
    netcdf.putVar(ncid, varido, oxygen);
    
    netcdf.putVar(ncid, varida, Talk);
    netcdf.putVar(ncid, variddic, DIC);
    
    %%close the netcdf file
    netcdf.close(ncid);
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____________________________________________________
% Modify Boundary file
% name of boundary file
bryname = [out_dir, bname];

%check to see if this file exists already, if not display warning and exit
if ~exist(bryname,'file')
    disp(['WARNING: Boundary file does not exist. Create a boundary before adding bio variables.']);
end

%% If the bry file does exist then proceed
if exist(bryname,'file')
    
    %open the netcdf file (assigns an id to it ato access it)
    fileattrib(bryname,'+w');
    ncid=netcdf.open(bryname, 'NC_WRITE');
    %% set the fill value to nothing so it will be zeros by default
    netcdf.setFill(ncid,'NC_NOFILL');
    
    %%  get id for variables of interest
    salt_time_id = netcdf.inqVarID(ncid,'salt_time');
    nsalt_id = netcdf.inqVarID(ncid,'salt_north');
    ssalt_id = netcdf.inqVarID(ncid,'salt_south');
    esalt_id = netcdf.inqVarID(ncid,'salt_east');
    wsalt_id = netcdf.inqVarID(ncid,'salt_west');
    ntemp_id = netcdf.inqVarID(ncid,'temp_north');
    stemp_id = netcdf.inqVarID(ncid,'temp_south');
    etemp_id = netcdf.inqVarID(ncid,'temp_east');
    wtemp_id = netcdf.inqVarID(ncid,'temp_west');
    %get the info about the salt including dimensions to copy to tracers
    [~, xtype, dimidsn, ~] = netcdf.inqVar(ncid,nsalt_id);
    [~, xtype, dimidss, ~] = netcdf.inqVar(ncid,ssalt_id);
    [~, xtype, dimidse, ~] = netcdf.inqVar(ncid,esalt_id);
    [~, xtype, dimidsw, ~] = netcdf.inqVar(ncid,wsalt_id);
    
    %Create time fields
    vartime = netcdf.getVar(ncid,salt_time_id,'double'); % basing bio time field off salt_time
    bry_nsalt = netcdf.getVar(ncid,nsalt_id,'double');
    bry_ssalt = netcdf.getVar(ncid,ssalt_id,'double');
    bry_esalt = netcdf.getVar(ncid,esalt_id,'double');
    bry_wsalt = netcdf.getVar(ncid,wsalt_id,'double');
    bry_ntemp = netcdf.getVar(ncid,ntemp_id,'double');
    bry_stemp = netcdf.getVar(ncid,stemp_id,'double');
    bry_etemp = netcdf.getVar(ncid,etemp_id,'double');
    bry_wtemp = netcdf.getVar(ncid,wtemp_id,'double');
    
    
    % Create NO3 field
    if strcmp(NO3_method,'PL_Salt')
        % Calculate NO3
        NO3n = make_NO3_field(NO3_method,bry_nsalt,bry_ntemp);
        NO3s = make_NO3_field(NO3_method,bry_ssalt,bry_stemp);
        NO3e = make_NO3_field(NO3_method,bry_esalt,bry_etemp);
        NO3w = make_NO3_field(NO3_method,bry_wsalt,bry_wtemp);
        % Calculate oxygen
        oxygenN = make_oxy_field(NO3_method,bry_nsalt);
        oxygenS = make_oxy_field(NO3_method,bry_ssalt);
        oxygenE = make_oxy_field(NO3_method,bry_esalt);
        oxygenW = make_oxy_field(NO3_method,bry_wsalt);
        [DIC_N, Talk_N]= make_DIC_field(NO3_method,bry_nsalt,bry_ntemp);
        [DIC_S, Talk_S] = make_DIC_field(NO3_method,bry_ssalt,bry_stemp);
        [DIC_E, Talk_E] = make_DIC_field(NO3_method,bry_esalt,bry_etemp);
        [DIC_W,Talk_W] = make_DIC_field(NO3_method,bry_wsalt,bry_wtemp);
        
        
    end
    
    %Create phytoplankton field
    phytoplanktonN = 0.01*ones(size(bry_nsalt)); % "Seed stock" of 0.01 microM of N
    phytoplanktonS = 0.01*ones(size(bry_ssalt)); % "Seed stock" of 0.01 microM of N
    phytoplanktonE = 0.01*ones(size(bry_esalt)); % "Seed stock" of 0.01 microM of N
    phytoplanktonW = 0.01*ones(size(bry_wsalt)); % "Seed stock" of 0.01 microM of N
    
    % Create zooplankton field
    zooplanktonN = 0.01*ones(size(bry_nsalt)); % "Seed stock" of 0.01 microM of N
    zooplanktonS = 0.01*ones(size(bry_ssalt)); % "Seed stock" of 0.01 microM of N
    zooplanktonE = 0.01*ones(size(bry_esalt)); % "Seed stock" of 0.01 microM of N
    zooplanktonW = 0.01*ones(size(bry_wsalt)); % "Seed stock" of 0.01 microM of N
    
    % Create detritus field
    detritusN = zeros(size(bry_nsalt)); % Initializing with no detritus
    detritusS = zeros(size(bry_ssalt)); % Initializing with no detritus
    detritusE = zeros(size(bry_esalt)); % Initializing with no detritus
    detritusW = zeros(size(bry_wsalt)); % Initializing with no detritus
    
    % Create detritus field
    LdetritusN = zeros(size(bry_nsalt)); % Initializing with no detritus
    LdetritusS = zeros(size(bry_ssalt)); % Initializing with no detritus
    LdetritusE = zeros(size(bry_esalt)); % Initializing with no detritus
    LdetritusW = zeros(size(bry_wsalt)); % Initializing with no detritus
    % Create PIC field
    PICN = zeros(size(bry_nsalt)); % Initializing with no detritus
    PICS = zeros(size(bry_ssalt)); % Initializing with no detritus
    PICE = zeros(size(bry_esalt)); % Initializing with no detritus
    PICW = zeros(size(bry_wsalt)); % Initializing with no detritus
    
    %%% WRITE VARIABLES TO NETCDF FILE %%%%%%%%%%%%%%%%%%%%%%
    %Define new time dimensions in clm file
    
    netcdf.reDef(ncid);
    dimidn=netcdf.defDim(ncid, 'NO3_time', length(vartime));
    dimidp=netcdf.defDim(ncid, 'phyt_time', length(vartime));
    dimidz=netcdf.defDim(ncid, 'zoop_time', length(vartime));
    dimidd=netcdf.defDim(ncid, 'detritus_time', length(vartime));
    dimidLd=netcdf.defDim(ncid, 'Ldetritus_time', length(vartime));
    dimidPIC=netcdf.defDim(ncid, 'CaCO3_time', length(vartime));
    dimido=netcdf.defDim(ncid, 'oxygen_time', length(vartime));
    dimida=netcdf.defDim(ncid, 'alk_time', length(vartime));
    dimiddic=netcdf.defDim(ncid, 'dic_time', length(vartime));
    
    
    %Define the variables based on salt dimensions
    varidnt=netcdf.defVar(ncid, 'NO3_time', 'double', dimidn);
    varidpt=netcdf.defVar(ncid, 'phyt_time', 'double', dimidp);
    varidzt=netcdf.defVar(ncid, 'zoop_time', 'double', dimidz);
    variddt=netcdf.defVar(ncid, 'detritus_time', 'double', dimidd);
    varidLdt=netcdf.defVar(ncid, 'Ldetritus_time', 'double', dimidLd);
    varidPICt=netcdf.defVar(ncid, 'CaCO3_time', 'double', dimidPIC);
    varidot=netcdf.defVar(ncid, 'oxygen_time', 'double', dimido);
    varidat=netcdf.defVar(ncid, 'alk_time', 'double', dimida);
    variddict=netcdf.defVar(ncid, 'dic_time', 'double', dimiddic);
    
    varidnN=netcdf.defVar(ncid, 'NO3_north',xtype,dimidsn);
    varidnS=netcdf.defVar(ncid, 'NO3_south',xtype,dimidss);
    varidnE=netcdf.defVar(ncid, 'NO3_east',xtype,dimidse);
    varidnW=netcdf.defVar(ncid, 'NO3_west',xtype,dimidsw);
    varidpN=netcdf.defVar(ncid, 'phytoplankton_north', xtype,dimidsn);
    varidpS=netcdf.defVar(ncid, 'phytoplankton_south', xtype,dimidss);
    varidpE=netcdf.defVar(ncid, 'phytoplankton_east', xtype,dimidse);
    varidpW=netcdf.defVar(ncid, 'phytoplankton_west', xtype,dimidsw);
    varidzN=netcdf.defVar(ncid, 'zooplankton_north',  xtype,dimidsn);
    varidzS=netcdf.defVar(ncid, 'zooplankton_south',  xtype,dimidss);
    varidzE=netcdf.defVar(ncid, 'zooplankton_east',  xtype,dimidse);
    varidzW=netcdf.defVar(ncid, 'zooplankton_west',  xtype,dimidsw);
    varidpN=netcdf.defVar(ncid, 'phyt_north', xtype,dimidsn);
    varidpS=netcdf.defVar(ncid, 'phyt_south', xtype,dimidss);
    varidpE=netcdf.defVar(ncid, 'phyt_east', xtype,dimidse);
    varidpW=netcdf.defVar(ncid, 'phyt_west', xtype,dimidsw);
    varidzN=netcdf.defVar(ncid, 'zoop_north',  xtype,dimidsn);
    varidzS=netcdf.defVar(ncid, 'zoop_south',  xtype,dimidss);
    varidzE=netcdf.defVar(ncid, 'zoop_east',  xtype,dimidse);
    varidzW=netcdf.defVar(ncid, 'zoop_west',  xtype,dimidsw);
    variddN=netcdf.defVar(ncid, 'detritus_north',  xtype,dimidsn);
    variddS=netcdf.defVar(ncid, 'detritus_south',  xtype,dimidss);
    variddE=netcdf.defVar(ncid, 'detritus_east',  xtype,dimidse);
    variddW=netcdf.defVar(ncid, 'detritus_west',  xtype,dimidsw);
    varidLdN=netcdf.defVar(ncid, 'Ldetritus_north',  xtype,dimidsn);
    varidLdS=netcdf.defVar(ncid, 'Ldetritus_south',  xtype,dimidss);
    varidLdE=netcdf.defVar(ncid, 'Ldetritus_east',  xtype,dimidse);
    varidLdW=netcdf.defVar(ncid, 'Ldetritus_west',  xtype,dimidsw);
    varidPICN=netcdf.defVar(ncid, 'CaCO3_north',  xtype,dimidsn);
    varidPICS=netcdf.defVar(ncid, 'CaCO3_south',  xtype,dimidss);
    varidPICE=netcdf.defVar(ncid, 'CaCO3_east',  xtype,dimidse);
    varidPICW=netcdf.defVar(ncid, 'CaCO3_west',  xtype,dimidsw);
    varidoN=netcdf.defVar(ncid, 'oxygen_north',  xtype,dimidsn);
    varidoS=netcdf.defVar(ncid, 'oxygen_south',  xtype,dimidss);
    varidoE=netcdf.defVar(ncid, 'oxygen_east',  xtype,dimidse);
    varidoW=netcdf.defVar(ncid, 'oxygen_west',  xtype,dimidsw);
    varidaN=netcdf.defVar(ncid, 'Talk_north',  xtype,dimidsn);
    varidaS=netcdf.defVar(ncid, 'Talk_south',  xtype,dimidss);
    varidaE=netcdf.defVar(ncid, 'Talk_east',  xtype,dimidse);
    varidaW=netcdf.defVar(ncid, 'Talk_west',  xtype,dimidsw);
    
    
    variddicN=netcdf.defVar(ncid, 'TIC_north',  xtype,dimidsn);
    variddicS=netcdf.defVar(ncid, 'TIC_south',  xtype,dimidss);
    variddicE=netcdf.defVar(ncid, 'TIC_east',  xtype,dimidse);
    variddicW=netcdf.defVar(ncid, 'TIC_west',  xtype,dimidsw);
    
    
    
    
    %add the desired attributes
    netcdf.putAtt(ncid,varidnt,'long_name','NO3 time');
    netcdf.putAtt(ncid,varidnt,'units','s');
    netcdf.putAtt(ncid,varidpt,'long_name','Phytoplankton time');
    netcdf.putAtt(ncid,varidpt,'units','s');
    netcdf.putAtt(ncid,varidzt,'long_name','Zooplankton time');
    netcdf.putAtt(ncid,varidzt,'units','s');
    netcdf.putAtt(ncid,variddt,'long_name','Detritus time');
    netcdf.putAtt(ncid,variddt,'units','s');
    netcdf.putAtt(ncid,varidLdt,'long_name','Large Detritus time');
    netcdf.putAtt(ncid,varidLdt,'units','s');
    netcdf.putAtt(ncid,varidPICt,'long_name','CaCO3 time');
    netcdf.putAtt(ncid,varidPICt,'units','s');
    netcdf.putAtt(ncid,varidot,'long_name','Oxygen time');
    netcdf.putAtt(ncid,varidot,'units','s');
    netcdf.putAtt(ncid,varidat,'long_name','alkalinity time');
    netcdf.putAtt(ncid,varidat,'units','s');
    netcdf.putAtt(ncid,variddict,'long_name','DIC time');
    netcdf.putAtt(ncid,variddict,'units','s');
    %
    %
    netcdf.putAtt(ncid,varidnN,'long_name','NO3 Northern BC');
    netcdf.putAtt(ncid,varidnN,'units','MicroMolar');
    netcdf.putAtt(ncid,varidnS,'long_name','NO3 Southern BC');
    netcdf.putAtt(ncid,varidnS,'units','MicroMolar');
    netcdf.putAtt(ncid,varidnE,'long_name','NO3 Eastern BC');
    netcdf.putAtt(ncid,varidnE,'units','MicroMolar');
    netcdf.putAtt(ncid,varidnW,'long_name','NO3 Western BC');
    netcdf.putAtt(ncid,varidnW,'units','MicroMolar');
    netcdf.putAtt(ncid,varidpN,'long_name','Phytoplankton Northern BC');
    netcdf.putAtt(ncid,varidpN,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidpS,'long_name','Phytoplankton Southern BC');
    netcdf.putAtt(ncid,varidpS,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidpE,'long_name','Phytoplankton Eastern BC');
    netcdf.putAtt(ncid,varidpE,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidpW,'long_name','Phytoplankton Western BC');
    netcdf.putAtt(ncid,varidpW,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidzN,'long_name','Zooplankton Northern BC');
    netcdf.putAtt(ncid,varidzN,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidzS,'long_name','Zooplankton Southern BC');
    netcdf.putAtt(ncid,varidzS,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidzE,'long_name','Zooplankton Eastern BC');
    netcdf.putAtt(ncid,varidzE,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidzW,'long_name','Zooplankton Western BC');
    netcdf.putAtt(ncid,varidzW,'units','MicroMolar N');
    netcdf.putAtt(ncid,variddN,'long_name','Detritus Northern BC');
    netcdf.putAtt(ncid,variddN,'units','MicroMolar N');
    netcdf.putAtt(ncid,variddS,'long_name','Detritus Southern BC');
    netcdf.putAtt(ncid,variddS,'units','MicroMolar N');
    netcdf.putAtt(ncid,variddE,'long_name','Detritus Eastern BC');
    netcdf.putAtt(ncid,variddE,'units','MicroMolar N');
    netcdf.putAtt(ncid,variddW,'long_name','Detritus Western BC');
    netcdf.putAtt(ncid,variddW,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidLdN,'long_name','Large Detritus Northern BC');
    netcdf.putAtt(ncid,varidLdN,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidLdS,'long_name','Large Detritus Southern BC');
    netcdf.putAtt(ncid,varidLdS,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidLdE,'long_name','Large Detritus Eastern BC');
    netcdf.putAtt(ncid,varidLdE,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidLdW,'long_name','Large Detritus Western BC');
    netcdf.putAtt(ncid,varidLdW,'units','MicroMolar N');
    netcdf.putAtt(ncid,varidPICN,'long_name','CaCO3 Northern BC');
    netcdf.putAtt(ncid,varidPICN,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidPICS,'long_name','CaCO3 Southern BC');
    netcdf.putAtt(ncid,varidPICS,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidPICE,'long_name','CaCO3 Eastern BC');
    netcdf.putAtt(ncid,varidPICE,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidPICW,'long_name','CaCO3 Western BC');
    netcdf.putAtt(ncid,varidPICW,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidoN,'long_name','Oxygen Northern BC');
    netcdf.putAtt(ncid,varidoN,'units','MicroMolar O');
    netcdf.putAtt(ncid,varidoS,'long_name','Oxygen Southern BC');
    netcdf.putAtt(ncid,varidoS,'units','MicroMolar O');
    netcdf.putAtt(ncid,varidoE,'long_name','Oxygen Eastern BC');
    netcdf.putAtt(ncid,varidoE,'units','MicroMolar O');
    netcdf.putAtt(ncid,varidoW,'long_name','Oxygen Western BC');
    netcdf.putAtt(ncid,varidoW,'units','MicroMolar O');
    netcdf.putAtt(ncid,varidaN,'long_name','Alk Northern BC');
    netcdf.putAtt(ncid,varidaN,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidaS,'long_name','Alk Southern BC');
    netcdf.putAtt(ncid,varidaS,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidaE,'long_name','Alk Eastern BC');
    netcdf.putAtt(ncid,varidaE,'units','MicroMolar C');
    netcdf.putAtt(ncid,varidaW,'long_name','Alk Western BC');
    netcdf.putAtt(ncid,varidaW,'units','MicroMolar C');
    
    netcdf.putAtt(ncid,variddicN,'long_name','TIC Northern BC');
    netcdf.putAtt(ncid,variddicN,'units','MicroMolar C');
    netcdf.putAtt(ncid,variddicS,'long_name','TIC Southern BC');
    netcdf.putAtt(ncid,variddicS,'units','MicroMolar O');
    netcdf.putAtt(ncid,variddicE,'long_name','TIC Eastern BC');
    netcdf.putAtt(ncid,variddicE,'units','MicroMolar C');
    netcdf.putAtt(ncid,variddicW,'long_name','TIC Western BC');
    netcdf.putAtt(ncid,variddicW,'units','MicroMolar C');
    
    %finish defining
    netcdf.endDef(ncid);
    
    %   Store the data
    netcdf.putVar(ncid,varidnt, vartime);
    netcdf.putVar(ncid,varidpt, vartime);
    netcdf.putVar(ncid,varidzt, vartime);
    netcdf.putVar(ncid,variddt, vartime);
    netcdf.putVar(ncid,varidLdt, vartime);
    netcdf.putVar(ncid,varidPICt, vartime);
    netcdf.putVar(ncid,varidot, vartime);
    netcdf.putVar(ncid,varidat, vartime);
    netcdf.putVar(ncid,variddict, vartime);
    
    netcdf.putVar(ncid, varidnN, NO3n);
    netcdf.putVar(ncid, varidnS, NO3s);
    netcdf.putVar(ncid, varidnE, NO3e);
    netcdf.putVar(ncid, varidnW, NO3w);
    netcdf.putVar(ncid, varidpN, phytoplanktonN);
    netcdf.putVar(ncid, varidpS, phytoplanktonS);
    netcdf.putVar(ncid, varidpE, phytoplanktonE);
    netcdf.putVar(ncid, varidpW, phytoplanktonW);
    netcdf.putVar(ncid, varidzN, zooplanktonN);
    netcdf.putVar(ncid, varidzS, zooplanktonS);
    netcdf.putVar(ncid, varidzE, zooplanktonE);
    netcdf.putVar(ncid, varidzW, zooplanktonW);
    netcdf.putVar(ncid, variddN, detritusN);
    netcdf.putVar(ncid, variddS, detritusS);
    netcdf.putVar(ncid, variddE, detritusE);
    netcdf.putVar(ncid, variddW, detritusW);
    netcdf.putVar(ncid, varidLdN, LdetritusN);
    netcdf.putVar(ncid, varidLdS, LdetritusS);
    netcdf.putVar(ncid, varidLdE, LdetritusE);
    netcdf.putVar(ncid, varidLdW, LdetritusW);
    netcdf.putVar(ncid, varidPICN, PICN);
    netcdf.putVar(ncid, varidPICS, PICS);
    netcdf.putVar(ncid, varidPICE, PICE);
    netcdf.putVar(ncid, varidPICW, PICW);
    netcdf.putVar(ncid, varidoN, oxygenN);
    netcdf.putVar(ncid, varidoS, oxygenS);
    netcdf.putVar(ncid, varidoE, oxygenE);
    netcdf.putVar(ncid, varidoW, oxygenW);
    
    netcdf.putVar(ncid, varidaN, Talk_N);
    netcdf.putVar(ncid, varidaS, Talk_S);
    netcdf.putVar(ncid, varidaE, Talk_E);
    netcdf.putVar(ncid, varidaW, Talk_W);
    
    netcdf.putVar(ncid, variddicN, DIC_N);
    netcdf.putVar(ncid, variddicS, DIC_S);
    netcdf.putVar(ncid, variddicE, DIC_E);
    netcdf.putVar(ncid, variddicW, DIC_W);
    
    
    %close the netcdf file
    netcdf.close(ncid);
    
    
end % END if exist(bryname)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %_____________________________________________________
% %Modify Initial conditions file
% %  if nargin>6 % If the name of an inifile is given, then add bio variables to the ini file
ininame = [out_dir, iname];

if exist(bryname,'file')
    
    
    %%open ini file to write to it
    fileattrib(ininame,'+w');
    ncidi=netcdf.open(ininame, 'NC_WRITE');
    ncidg=netcdf.open(gridname, 'NC_NOWRITE');
    %     %% set the fill value to nothing so it will be zeros by default
    netcdf.setFill(ncidi,'NC_NOFILL');
    %
    %     %%get id for variables of interest
    salt_time_id = netcdf.inqVarID(ncidi,'ocean_time');
    salt_id = netcdf.inqVarID(ncidi,'salt');
    temp_id = netcdf.inqVarID(ncidi,'temp');
    %     zeta_id = netcdf.inqVarID(ncidi,'zeta');
    %     %h_id = netcdf.inqVarID(ncidg,'h');
    %
    %     %%get the info about the salt including dimensions to copy to tracers
    [~, xtype, dimids, ~] = netcdf.inqVar(ncidi,salt_id);
    %
    %     %%Create time fields
    vartime = netcdf.getVar(ncidi,salt_time_id,'double'); % basing bio time field off salt_time
    ini_salt = netcdf.getVar(ncidi,salt_id,'double');
    ini_temp = netcdf.getVar(ncidi,temp_id,'double');
    %ini_zeta = netcdf.getVar(ncidi,zeta_id,'double');
    %     %h_grid=netcdf.getVar(ncidg,h_id,'double');
    %
    %     %%Calc z to pass to NO3 and O2 generating schemes
    %     %ini_dep=roms_z(h_grid,ini_zeta);
    %
    %     %%Create NO3 field
    if strcmp(NO3_method,'PL_Salt')
        %Calculate NO3
        NO3 = make_NO3_field(NO3_method,ini_salt);
        % Now set nitrate in Salish Sea
        % Read in the grid matrices
        lon_rho = nc_varget(gname,'lon_rho')';
        lat_rho = nc_varget(gname,'lat_rho')';
        [iB1,jB1] = find(and(lat_rho>=49,lon_rho>=-125.25)); %Box 1
        [iB2,jB2] = find(and(and(lat_rho<49,lat_rho>=48.5),lon_rho>=-124)); %Box 2
        [iB3,jB3] = find(and(and(lat_rho<48.5,lat_rho>=47),lon_rho>=-123.5)); %Box 3
        salt_corr = 0.9276*ini_salt + 2.746; % NCOM salinity correction
        NO3(iB1,jB1,:) = 3.26*salt_corr(iB1,jB1,:) -76.44; % line from Neil's fit to Eastern SJDF data
        NO3(iB2,jB2,:) = 3.26*salt_corr(iB2,jB2,:) -76.44; % line from Neil's fit to Eastern SJDF data
        NO3(iB3,jB3,:) = 3.26*salt_corr(iB3,jB3,:) -76.44; % line from Neil's fit to Eastern SJDF data
        clear index
        
        index = find(NO3 < 0);
        NO3(index) = 0;
        clear index
        
        %Calculate oxygen
        oxygen = make_oxy_field(NO3_method,ini_salt);
        [DIC, Talk] = make_DIC_field(NO3_method,ini_salt,ini_temp);
        
    end
    
    % Create phytoplankton field
    phytoplankton = 0.01*ones(size(ini_salt)); % "Seed stock" of 0.01 microM of N
    
    % Create zooplankton field
    zooplankton = 0.01*ones(size(ini_salt)); % "Seed stock" of 0.01 microM of N
    
    % Create detritus field
    detritus = zeros(size(ini_salt)); % Initializing with no detritus
    
    % Create detritus field
    Ldetritus = zeros(size(ini_salt)); % Initializing with no detritus
    PIC = zeros(size(ini_salt)); % Initializing with no PIC
    
    %add dimensions to ini file
    % turn on defining
    netcdf.reDef(ncidi);
    
    
    %define the variable based on the salt dimensions and type
    varidn = netcdf.defVar(ncidi,'NO3',xtype,dimids);
    netcdf.putAtt(ncidi,varidn,'long_name','NO3');
    netcdf.putAtt(ncidi,varidn,'units','MicroMolar');
    
    varidp = netcdf.defVar(ncidi,'phytoplankton',xtype,dimids);
    netcdf.putAtt(ncidi,varidp,'long_name','phytoplankton');
    netcdf.putAtt(ncidi,varidp,'units','MicroMolar N');
    
    varidz = netcdf.defVar(ncidi,'zooplankton',xtype,dimids);
    netcdf.putAtt(ncidi,varidz,'long_name','zooplankton');
    netcdf.putAtt(ncidi,varidz,'units','MicroMolar N');
    
    varidd = netcdf.defVar(ncidi,'detritus',xtype,dimids);
    netcdf.putAtt(ncidi,varidd,'long_name','detritus');
    netcdf.putAtt(ncidi,varidd,'units','MicroMolar N');
    
    varidLd = netcdf.defVar(ncidi,'Ldetritus',xtype,dimids);
    netcdf.putAtt(ncidi,varidLd,'long_name','Ldetritus');
    netcdf.putAtt(ncidi,varidLd,'units','MicroMolar N');
    
    varidPIC = netcdf.defVar(ncidi,'CaCO3',xtype,dimids);
    netcdf.putAtt(ncidi,varidPIC,'long_name','CaCO3');
    netcdf.putAtt(ncidi,varidPIC,'units','MicroMolar C');
    
    varido = netcdf.defVar(ncidi,'oxygen',xtype,dimids);
    netcdf.putAtt(ncidi,varido,'long_name','oxygen');
    netcdf.putAtt(ncidi,varido,'units','MicroMolar');
    
    varida = netcdf.defVar(ncidi,'alkalinity',xtype,dimids);
    netcdf.putAtt(ncidi,varida,'long_name','alkalinity');
    netcdf.putAtt(ncidi,varida,'units','MicroMolar');
    
    variddic = netcdf.defVar(ncidi,'TIC',xtype,dimids);
    netcdf.putAtt(ncidi,variddic,'long_name','DIC');
    netcdf.putAtt(ncidi,variddic,'units','MicroMolar');
    
    %finish defining
    netcdf.endDef(ncidi);
    
    
    %put the fields into the ini file
    %Store the data
    netcdf.putVar(ncidi, varidn, NO3);
    netcdf.putVar(ncidi, varidp, phytoplankton);
    netcdf.putVar(ncidi, varidz, zooplankton);
    netcdf.putVar(ncidi, varidd, detritus);
    netcdf.putVar(ncidi, varidLd, Ldetritus);
    netcdf.putVar(ncidi, varidPIC, PIC);
    netcdf.putVar(ncidi, varido, oxygen);
    netcdf.putVar(ncidi, varida, Talk);
    netcdf.putVar(ncidi, variddic, DIC);
    
    %close the netcdf file
    netcdf.close(ncidi);
    
end % INI if statement

% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Modifying rivers.nc input file
% % % % name of river forcing file
rivname = [out_dir, rname];
% % % % If the rivers.nc file does exist then proceed
if exist(rivname,'file')
    % % % %     % Create time fields
    vartime = nc_varget(rivname,'river_time'); % basing bio time field off river_time
    river_salt = nc_varget(rivname,'river_salt');
    %     % Create NO3 field
    %     % Create NO3 field
    %     % Load file with river nutrient climatology (Columbia and Fraser only)
    load('/home/disk/edward/siedlesa/River_TotalN_clm.mat')
    %     % Add one day to climatology file for 2005 because Sarah's new rivers
    %     % file has 366 days:
    CR_TotalN_clm = [CR_TotalN_clm; CR_TotalN_clm(end)];
    FR_TotalN_clm = [FR_TotalN_clm; FR_TotalN_clm(end)];
    %
    river_NO3 = 5*ones(size(river_salt)); % initialize nitrate array with 5 microMolar N
    %     % River names (ID#): (1) Skagit, (2) Snohomish, (3) Stilliguamish, (4)
    %     % Columbia - 2 inputs, (5) Puyallup, (6) Duwamish, (7) Nisqually,
    %     % (8)Deschutes, (9) Skagit-South, (10) Skokomish, (11) Duckabush, (12)
    %     % Fraser, (13) Dosewallips, (14) Hammahamma, (15) Cedar, (16) Nooksack,
    %     % (17) Samish
    %     % Set Columbia and Fraser River nutrients using seasonal climatology derived
    %     % from data
    river_NO3(:,:,4) = repmat(CR_TotalN_clm,1,size(river_salt,2)); % Columbia input 1
    river_NO3(:,:,5) = repmat(CR_TotalN_clm,1,size(river_salt,2)); % Columbia input 2
    river_NO3(:,:,12) = repmat(FR_TotalN_clm,1,size(river_salt,2)); % Fraser
    %
    %     river_NO3 = 5*ones(size(river_salt)); % setting the nitrate level in all rivers to 5 microMolar
    % %
    % %     % Create phytoplankton field
    river_phytoplankton = 0.01*ones(size(river_salt)); % "Seed stock" of 0.01 microM of N
    % %
    % %     % Create zooplankton field
    river_zooplankton = 0.01*ones(size(river_salt)); % "Seed stock" of 0.01 microM of N
    %
    river_detritus = zeros(size(river_salt)); % Initializing with a small amount of detritus = 1 microMolar
    %
    river_Ldetritus = zeros(size(river_salt)); % Initializing with a small amount of detritus = 1 microMolar
    river_CaCO3 = zeros(size(river_salt)); % Initializing with a small amount of detritus = 1 microMolar
    %
    %      % Create oxygen field
    river_oxygen = 350*ones(size(river_salt)); % setting the oxygen level in all rivers to 350 microMolar
    river_TIC = 1700*ones(size(river_salt)); % setting the DIC level
    river_alkalinity = 900*ones(size(river_salt)); % setting the alkalinity
    
    %
    %     %%%% WRITE VARIABLES TO NETCDF FILE %%%%%%%%%%%%%%%%%%%%%%
    %
    %     % Add the variables
    varstruct.Name = 'river_NO3'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff NO3','MicroMolar'});
    nc_addvar(rivname, varstruct);
    %
    varstruct.Name = 'river_phytoplankton'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff phytoplankton concentration','MicroMolar N'});
    nc_addvar(rivname, varstruct);
    %
    varstruct.Name = 'river_zooplankton'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff zooplankton concentration','MicroMolar N'});
    nc_addvar(rivname, varstruct);
    %
    varstruct.Name = 'river_detritus'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff detritus concentration','MicroMolar N'});
    nc_addvar(rivname, varstruct);
    %
    varstruct.Name = 'river_Ldetritus'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff Ldetritus concentration','MicroMolar N'});
    nc_addvar(rivname, varstruct);
    varstruct.Name = 'river_CaCO3'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff CaCO3 concentration','MicroMolar C'});
    nc_addvar(rivname, varstruct);
    %
    varstruct.Name = 'river_oxygen'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff oxygen concentration','millimole/m3'});
    nc_addvar(rivname, varstruct);
    varstruct.Name = 'river_TIC'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff TIC concentration','millimole/m3'});
    nc_addvar(rivname, varstruct);
    varstruct.Name = 'river_alkalinity'; varstruct.Dimension = {'river_time','s_rho','river'};
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{'river runoff alkalinity concentration','millimole/m3'});
    nc_addvar(rivname, varstruct);
    %
    %     % Store the data
    nc_varput(rivname, 'river_NO3', river_NO3);
    nc_varput(rivname, 'river_phytoplankton', river_phytoplankton);
    nc_varput(rivname, 'river_zooplankton', river_zooplankton);
    nc_varput(rivname, 'river_detritus', river_detritus);
    nc_varput(rivname, 'river_Ldetritus', river_Ldetritus);
    nc_varput(rivname, 'river_CaCO3', river_CaCO3);
    nc_varput(rivname, 'river_oxygen', river_oxygen);
    nc_varput(rivname, 'river_TIC', river_TIC);
    nc_varput(rivname, 'river_alkalinity', river_alkalinity);
end
%

disp('DONE')

toc

