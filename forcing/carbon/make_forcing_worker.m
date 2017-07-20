function make_forcing_worker(gridname, tag, date_string, run_type, indir, h0, h1, outdir)

%% make_forcing_worker.m
%
% ****************** for carbon ***********************
%
% This adds pH and Aragonite saturation state to a ROMS history file.

addpath('../../alpha'); Ldir = Lstart(gridname, tag);
start_time = datenum(now);

%% carbon-specific code

addpath('../../shared/seawater');
addpath('../../shared'); % gives access to CO2SYS_PM.m

%% create the output time vector

% if strcmp(run_type,'backfill')
%     hr_vec = 1:24;
%     %hr_vec = 1:2;
%     % these correspond to the hour of the day, but note that the his file
%     % number is one greater
% elseif strcmp(run_type,'forecast')
%     hr_vec = 1:72;
%     % this gives 72 hourly values
%     % (three days, including endpoint, but not hour zero)
% elseif strcmp(run_type,'low_passed')
%     hr_vec = 1;
%     % a new option, that just works on low_passed.nc
% end

h0 = str2num(h0); h1 = str2num(h1);
hr_vec = [h0:h1];

%% loop over the hours

for tt = 1:length(hr_vec)
    
    if strcmp(run_type,'low_passed')
        fn = 'low_passed.nc';
    else
        hr_num = hr_vec(tt);
        hr_string = ['0000', num2str(hr_num)];
        hr_string = hr_string(end-3:end);
        fn = ['ocean_his_', hr_string, '.nc'];
    end
    
    infile = [indir,fn];
    disp(' ')
    disp(infile)
    
    % driver code for CO2SYS
    % we make use mainly of the nc... convenience functions that are
    % now part of MATLAB
    
    % note that these are packed (x, y, z, t)
    TIC = ncread(infile, 'TIC');
    alkalinity= ncread(infile, 'alkalinity');
    salt = ncread(infile, 'salt');
    temp = ncread(infile, 'temp');
    
    % other stuff we need
    %
    % these are packed (t, z, y, x) or some appropriate subset of that
    [S] = Z_get_S(infile);
    h = nc_varget(infile,'h');
    [z_rho] = Z_s2z_rho(h, 0*h, S);
    % so we permute this to make it like salt, etc.
    z_rho = permute(z_rho, [3, 2, 1]);
    
    % % prepare vectors for calculation (no need to strip NaNs)
    DIC0 = TIC(:);
    ALK0 = alkalinity(:);
    SALT = salt(:);
    THETA = temp(:); % potential temperature
    Z_RHO = z_rho(:);
    PRES = sw_pres(-Z_RHO, 45); % decibars, not Pa
    TEMP = sw_ptmp(SALT, THETA, 0, PRES); % in situ temperature
    DEN = sw_dens(SALT, TEMP, PRES); % in situ density
    % Convert from umol/m3 to umol/kg?
    % Are we really supposed to be using in situ density for this?
    % Recall that ROMS only considers density to be compressible
    % for the calculation of pressure gradients, not for volume
    % conservation
    DIC = 1000*DIC0./DEN;
    ALK = 1000*ALK0./DEN;
    
    % **** FIX BAD VALUES ****
    % 11/28/2016 there were two near-zero values of TIC along the
    % north row, and these caused CO2SYS to choke.
    % The 100 value is arbitrary, but seemed reasonable to me.
    DIC(DIC<100) = NaN;
    ALK(ALK<100) = NaN;
    
    % prepare to call CO2SYS
    
    % We break up the calculation into chunks "nn" long, to avoid
    % overloadiang the function call (performance degrades as nn -> NN).
    %
    % RESULT: the CO2SYS calculation took about 10 seconds on my mac
    % for a single cascadia1_base_lobio1 history file, with nn = 100000.
    nn = 100000;
    
    % pre-allocate
    NN = length(SALT);   
    PH = NaN * SALT;
    ARAG = NaN * SALT;
    
    tic
    
    % for debugging use nn = 2 and i0 = 716241
    for i0 = 1:nn:NN
        
        %disp(['Item ',num2str(i0),' out of ',num2str(NN)])
        
        i1 = i0 + nn - 1;
        if i1 > NN; i1 = NN; end;
        
        % Assume IN means in situ, and OUT means at the surface (changes pressure)
        PAR1 = ALK(i0:i1); %  (some unit) : scalar or vector of size n
        PAR2 = DIC(i0:i1); %  (some unit) : scalar or vector of size n
        PAR1TYPE = 1; %      () : scalar or vector of size n (*)
        PAR2TYPE = 2; %     () : scalar or vector of size n (*)
        SAL = SALT(i0:i1); %            () : scalar or vector of size n
        TEMPIN = TEMP(i0:i1); %  (degr. C) : scalar or vector of size n
        TEMPOUT = THETA(i0:i1); % (degr. C) : scalar or vector of size n
        PRESIN = PRES(i0:i1); %     (dbar) : scalar or vector of size n
        PRESOUT = 0; %   (dbar) : scalar or vector of size n
        SI = 50; %    (umol/kgSW) : scalar or vector of size n
        PO4 = 2; %   (umol/kgSW) : scalar or vector of size n
        pHSCALEIN = 1; %        : scalar or vector of size n (**)
        K1K2CONSTANTS = 10; %     : scalar or vector of size n (***)
        KSO4CONSTANTS = 1; %    : scalar or vector of size n (****)
        
        % Do the calculation. See CO2SYS's help for syntax and output format
        A=CO2SYS_PM(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
            SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,...
            K1K2CONSTANTS,KSO4CONSTANTS);
        PH(i0:i1) = A(:,3);
        ARAG(i0:i1) = A(:,16);
        
    end
    
    dt = toc;
    disp(['Took ',num2str(round(dt)),' seconds'])
    
    % reshape to be like the history file fields
    PH = reshape(PH, [size(salt),1]);
    ARAG = reshape(ARAG, [size(salt),1]);
    
    % Write to the NetCDf history files
    
    outvar_list = {'PH', 'ARAG'};
    
    for ii = 1:length(outvar_list)
        vn = outvar_list{ii};
        eval(['this_var = ',vn,';']);
        try
            % variable already exists
            vn_info = ncinfo(infile,vn);
            disp(['filling existing ',vn])
            ncwrite(infile, vn, this_var);
        catch exception
            % variable does not exist yet
            if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
                disp(['creating new ',vn])
                sch =  ncinfo(infile, 'salt');
                sch.Name = vn;
                % Ae should also overwrite the long_name, but I am
                % not going to bother.  It is in
                % sch.Attribute(1).Value.
                ncwriteschema(infile, sch);
                ncwrite(infile, vn, this_var);
            end
        end
        
    end
    
end % end of hour loop

%% things to use for checking result

% Actually I'm not sure how we check the result for success.  We don't 
% make new files, just change existing ones.

%% Final output
datestr_format = 'yyyy.mm.dd HH:MM:SS';
end_time = datenum(now);
fid = fopen([outdir,'Info/process_status.csv'],'w');
fprintf(fid,'%s\n',['start_time,',datestr(start_time, datestr_format)]);
fprintf(fid,'%s\n',['end_time,',datestr(end_time, datestr_format)]);
% for now assume it succeeded.
fprintf(fid,'%s\n','result,success');
fclose(fid);

