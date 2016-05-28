function run_co2sys(indir, moor_file, in3)
% indir = '/Users/PM5/Documents/LiveOcean_output/moor/'
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

addpath('../shared/seawater')

ncid = netcdf.open([indir,moor_file], 'WRITE')
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

var2get_list = {'TIC','alkalinity','salt','temp'};
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

s_rho_id = netcdf.inqDimID(ncid,'s_rho')
ot_id = netcdf.inqDimID(ncid,'ocean_time')

%netcdf.close(ncid)


if 1
% % prepare vectors for calculation (no need to strip NaNs)
Dic = TIC(:);
Ta = alkalinity(:);
Salt = salt(:);
Temp = temp(:);
Z_rho = z_rho(:);
Lat_rho = lat_rho(:);
Pres = sw_pres(-Z_rho,lat_rho);


%% carbon calc

% We break up the calculation into chunks "nn" long, to avoid
% overwhelming the function call (performance degrades as nn -> NN).

NN = length(Salt);
nn = 1000;

PH = NaN * Salt;
ARAG = NaN * Salt;

tic

for i0 = 1:nn:NN
    
    i1 = i0 + nn - 1;
    if i1 > NN; i1 = NN; end;
    
    disp(['i0 = ',num2str(i0),' i1 = ',num2str(i1)])
    
    par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
    par1     = Ta(i0:i1); % value of the first parameter
    par2type =    2; % The first parameter supplied is of type "1", which is "DIC"
    par2     = Dic(i0:i1); % value of the second parameter, which is a long vector of different DIC's!
    sal      =  Salt(i0:i1); % Salinity of the sample
    tempin   =   Temp(i0:i1); % Temperature at input conditions
    presin   =    Pres(i0:i1); % Pressure    at input conditions
    tempout  =    0; % Temperature at output conditions - doesn't matter in this example
    presout  =    0; % Pressure    at output conditions - doesn't matter in this example
    sil      =   50; % Concentration of silicate  in the sample (in umol/kg)
    po4      =    2; % Concentration of phosphate in the sample (in umol/kg)
    pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
    k1k2c    =    10; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    
    % Do the calculation. See CO2SYS's help for syntax and output format
    A = CO2SYS_PM(par1,par2,par1type,par2type,sal, ...
        tempin,tempout,presin,presout, ...
        sil,po4,pHscale,k1k2c,kso4c);
    PH(i0:i1) = A(:,18);
    ARAG(i0:i1) = A(:,16);
    clear A par1 par2 sal tempin presin
    
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
        varid = netcdf.defVar(ncid,var2put,'double',[ot_id, s_rho_id]);
        netcdf.putVar(ncid,varid,var_temp)
    end
end

netcdf.close(ncid);

end




