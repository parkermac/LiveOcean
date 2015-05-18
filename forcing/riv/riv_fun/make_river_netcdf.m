function make_river_netcdf(MP, LP, S, rout, riv_file_out, river_time, ...
    Qr_flow, T_riv, uid, initialize)

% called by make_forcing_worker.m

%% ------------------------------------------------------------------------
% Create a netcdf file that contains river forcing data for ROMS.
% Forcing data consists of:
% 'river_Xposition'  -   'river runoff  XI-positions at RHO-points'
% 'river_Eposition'  -   'river runoff ETA-positions at RHO-points'
% 'river_direction'  -   'river runoff direction'
% 'river_Vshape'     -   'river runoff mass transport vertical profile'
% 'river_transport'  -   'river runoff mass transport'
% 'river_flag'       -   'river runoff flag'
% 'river_temp'       -   'river runoff potential temperature'
% 'river_salt'       -   'river runoff salinity'
%
% In cppdefs.h you should have
% #define TS_PSOURCE
% #define UV_PSOURCE
% #undef  ANA_PSOURCE
% #undef  RIVER_SEDIMENT

%********************************************************
% calc some grid stuff here
%********************************************************
Lm = LP-2;
Mm = MP-2;
N = S.N;
L  = Lm+1;
M  = Mm+1;

% number of times
num_river_times = length(river_time);

% now initialize river locations, directions, and flags
river_Xposition = rout(1).X(:);
river_Eposition = rout(1).Y(:);
river_direction = rout(1).D(:);
sign = rout(1).sign(:);
numcells = length(rout(1).X)*ones(length(rout(1).X),1);
river_ID = ones(length(rout(1).X),1);
numr = length(uid); %
if(numr > 1) %more than 1 river
    for j=2:numr
        river_Xposition = [river_Xposition;rout(j).X(:)];
        river_Eposition = [river_Eposition;rout(j).Y(:)];
        river_direction = [river_direction;rout(j).D(:)];
        sign = [sign;rout(j).sign(:)];
        numcells = [numcells; length(rout(j).X)*ones(length(rout(j).X),1)];
        river_ID = [river_ID; j*ones(length(rout(j).X),1)];
    end
end
Nsources = length(river_Xposition); % could be different than # of rivers
river_flag = 3*ones(Nsources,1);
% "3" specifies temp and sal at each point source

% initialize river shape.
for i=1:Nsources
    river_Vshape(:,i) = flipud(linspace(2/N,0,N)');
    %linear decay from surface value to 0, fixed by sng 7/2011
end

% ramp river flow if needed
river_transport = nan(num_river_times,Nsources);
for i=1:Nsources
    fac = ones(size(river_time));
    if initialize==1
        fac(river_time<43200) = 1.0+tanh((river_time(time)-43200.0)/43200.0);
        disp('ramping up river flow');
    end
    Qi = Qr_flow(:,river_ID(i));
    Q = Qi./numcells(i); %divide this river Q by number of cells
    river_transport(:,i)=sign(i).*fac.*Q; %make negative for S or W
end

river_temp = nan(num_river_times,N,Nsources);
river_salt = zeros(num_river_times,N,Nsources);
for k=1:N
    for i=1:Nsources
        river_temp(:,k,i) = T_riv(:,river_ID(i));
    end
end

NAT = 2;    % assume temp + salt are active
NT = NAT;   % total number of tracers.

%%
my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
disp(['*** creating ',riv_file_out,' ***']);
nc_create_empty ( riv_file_out, my_mode );
nc_padheader (riv_file_out , 30000 );

% Global attributes:
disp(' ## Defining Global Attributes...')
nc_attput(riv_file_out,nc_global,'history', ...
    ['Created by "' mfilename '" on ' datestr(now)]);
nc_attput(riv_file_out,nc_global,'type', ...
    'Initialization file from create_roms_init.m');

% Dimensions:
disp(' ## Defining Dimensions...')

nc_add_dimension(riv_file_out, 'xi_rho', LP);
nc_add_dimension(riv_file_out, 'eta_rho', MP);
nc_add_dimension(riv_file_out, 'xi_psi', L);
nc_add_dimension(riv_file_out, 'eta_psi', M);
nc_add_dimension(riv_file_out, 'xi_u', L);
nc_add_dimension(riv_file_out, 'eta_u', MP);
nc_add_dimension(riv_file_out, 'xi_v', LP);
nc_add_dimension(riv_file_out, 'eta_v', M);

nc_add_dimension(riv_file_out, 's_rho', N);
nc_add_dimension(riv_file_out, 'tracer', NT);
nc_add_dimension(riv_file_out, 's_w', N+1);

nc_add_dimension(riv_file_out, 'one', 1);
nc_add_dimension(riv_file_out, 'two', 2);
nc_add_dimension(riv_file_out, 'river', Nsources);
nc_add_dimension(riv_file_out, 'river_time', num_river_times);

% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

varstruct.Name = 'theta_b';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate bottom control parameter'];
units = '1';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'theta_s';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate surface control parameter'];
units = '1';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Tcline';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate surface/bottom layer width'];
units = 'meter';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'hc';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate parameter, critical depth'];
units = 'meter';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Cs_r';
varstruct.Dimension = {'s_rho'};
long_name = ['S-coordinate stretching curves at RHO-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'Cs_r, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Cs_w';
varstruct.Dimension = {'s_w'};
long_name = ['S-coordinate stretching curves at W-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'Cs_w, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'sc_r';
varstruct.Dimension = {'s_rho'};
long_name = ['S-coordinate at RHO-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'sc_r, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'sc_w';
varstruct.Dimension = {'s_w'};
long_name = ['S-coordinate at W-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'sc_w, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river';
varstruct.Dimension = {'river'};
long_name = ['river_runoff identification number'];
units = 'nondimensional';
field = 'num_rivers, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_time';
varstruct.Dimension = {'river_time'};
long_name = ['river_time'];
units = 'seconds';
field = 'river_time, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Xposition';
varstruct.Dimension = {'river'};
long_name = ['river runoff  XI-positions at RHO-points'];
units = 'scalar';
field = 'river runoff XI position, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Eposition';
varstruct.Dimension = {'river'};
long_name = ['river runoff  ETA-positions at RHO-points'];
units = 'scalar';
field = 'river runoff ETA position, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_direction';
varstruct.Dimension = {'river'};
long_name = ['river runoff direction, XI=0, ETA>0'];
units = 'scalar';
field = 'river runoff direction, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Vshape';
varstruct.Dimension = {'s_rho','river'};
long_name = ['river runoff mass transport vertical profile'];
units = 'scalar';
field = 'river runoff vertical profile, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_transport';
varstruct.Dimension = {'river_time','river'};
long_name = ['river runoff mass transport'];
units = 'meter^3/s';
field = 'river runoff mass transport, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_flag';
varstruct.Dimension = {'river'};
long_name = ['river flag, 1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio'];
units = 'nondimensional';
field = 'river flag, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_temp';
varstruct.Dimension = {'river_time','s_rho','river'};
long_name = ['river runoff potential temperature'];
units = 'Celsius';
field = 'river temp, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_salt';
varstruct.Dimension = {'river_time','s_rho','river'};
long_name = ['river runoff salinity'];
units = 'psu';
field = 'river salinity, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')
nc_varput(riv_file_out, 'theta_s', S.theta_s);
nc_varput(riv_file_out, 'theta_b', S.theta_b);
nc_varput(riv_file_out, 'Tcline', S.tcline);
nc_varput(riv_file_out, 'Cs_r', S.Cs_r);
nc_varput(riv_file_out, 'Cs_w', S.Cs_w);
nc_varput(riv_file_out, 'sc_w', S.s_w);
nc_varput(riv_file_out, 'sc_r', S.s_rho);
nc_varput(riv_file_out, 'hc', S.hc);
nc_varput(riv_file_out, 'river', river_ID);

nc_varput(riv_file_out, 'river_time', river_time);
nc_varput(riv_file_out, 'river_Xposition', river_Xposition);
nc_varput(riv_file_out, 'river_Eposition', river_Eposition);
nc_varput(riv_file_out, 'river_direction', river_direction);
nc_varput(riv_file_out, 'river_Vshape', river_Vshape);
nc_varput(riv_file_out, 'river_transport', river_transport);
nc_varput(riv_file_out, 'river_flag', river_flag);
nc_varput(riv_file_out, 'river_temp', river_temp);
nc_varput(riv_file_out, 'river_salt', river_salt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
