function run_co2sys(indir, moor_file, in3)
% indir = '~/Documents/LiveOcean_output/moor/'
% moor_file = 'cascadia1_base_lobio1_NH10_low_pass_2013.06.02_2014.12.30.nc'

% test of co2sys with ROMS fields
%
% RESULT: the code would hang using CO2SYS.m, and the problem appeared
% to be the iterations performed in the function CalculatepHfromTATC().
% This has an open-ended while loop that was causing the problem. By
% increasing the value of pHTol from 0.0001 to 0.001 the performance
% improved hugely.
%
% This change of pHTol is implemented in CO2SYS_PM.m, along with
% a few lines to print the number of iterations.
% All edits can be found by searching for "MacCready".
%
% performance notes:
%
% takes 4 sec with nn = 1000 and pHTol = 0.001
% for a full surface field

% *** Need to rewrite following the forcing/carbon version which is
% cleaner. ***

addpath('../shared/seawater')
addpath('../shared'); % gives access to CO2SYS_PM.m

ncid = netcdf.open([indir,moor_file], 'WRITE');
% get the full listing of variables and their info
varid_list = netcdf.inqVarIDs(ncid);
varname_list = cell(1,length(varid_list));
dimids_list = cell(1,length(varid_list));
for ii = 1:length(varid_list)
    varid = varid_list(ii);
    [varname_list{ii},~,dimids_list{ii},~] = ...
        netcdf.inqVar(ncid,varid);
end
% for ii = 1:length(varname_list)
%     disp(varname_list{ii})
% end

var2get_list = {'TIC','alkalinity','salt','temp','z_rho'};
for vv = 1:length(var2get_list)
    var2get = var2get_list{vv};
    tf = strcmp(var2get,varname_list);
    if sum(tf)==1
        % do this if the variable exists
        index = find(tf);
        varname = varname_list{index};
        varid = varid_list(index);
        var = netcdf.getVar(ncid,varid,'double');
        % NOTE that using netcdf.getVar the fields are packed (x,y,z,t) so
        % the lat and lon matrices are packed (x,y), but this is confusing
        % so we will convert EVERYTHING back to (y,x)
        var = var'; % now stored as var(y,x)
        eval([varname,' = var;']);
    else
        % do this if the variable is absent
        disp(['** ',var2get,' not found!']);
    end
end

s_rho_id = netcdf.inqDimID(ncid,'s_rho');
ot_id = netcdf.inqDimID(ncid,'ocean_time');

if 1
    % % prepare vectors for calculation (no need to strip NaNs)
    DIC0 = TIC(:);
    ALK0 = alkalinity(:);
    SALT = salt(:);
    THETA = temp(:);
    Z_RHO = z_rho(:);
    PRES = sw_pres(-Z_RHO, 45); % decibars, not Pa
    TEMP = sw_ptmp(SALT, THETA, 0, PRES);
    DEN = sw_dens(SALT, TEMP, PRES);
    DIC = 1000*DIC0./DEN;
    ALK = 1000*ALK0./DEN;
        
    %% carbon calc
    
    % We break up the calculation into chunks "nn" long, to avoid
    % overwhelming the function call (performance degrades as nn -> NN).
    
    NN = length(SALT);
    nn = 1000;
    
    PH = NaN * SALT;
    ARAG = NaN * SALT;
    
    tic
    
    for i0 = 1:nn:NN
        
        i1 = i0 + nn - 1;
        if i1 > NN; i1 = NN; end;
                
        % Assume in mean in situ, and out means at the surface (changes pressure)
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
        clear A PAR1 PAR2 SAL TEMPIN PRESIN
        
    end
    
    dt = toc;
    disp(['Took ',num2str(round(dt)),' seconds'])
    
    %% write output
    PH = reshape(PH, size(salt))';
    ARAG = reshape(ARAG, size(salt))';
    
    var2put_list = {'PH','ARAG'};
    for vv = 1:length(var2put_list)
        var2put = var2put_list{vv};
        eval(['var_temp = ', var2put,';']);
        
        tf = strcmp(var2put,varname_list);
        if sum(tf)==1
            % do this if the variable exists
            disp(['** ',var2put,' found']);
            index = find(tf);
            %varname = varname_list{index};
            varid = varid_list(index);
            netcdf.putVar(ncid,varid,var_temp)
        else
            % do this if the variable is absent
            disp(['** ',var2put,' not found']);
            varid = netcdf.defVar(ncid,var2put,'double',[s_rho_id, ot_id]);
            netcdf.putVar(ncid,varid,var_temp)
        end
    end
end

netcdf.close(ncid);





