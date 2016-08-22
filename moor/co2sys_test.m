% code to test CO2SYS

par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     = 2325; % value of the first parameter
par2type =    2; % The first parameter supplied is of type "1", which is "DIC"
par2     = 2400; % value of the second parameter, which is a long vector of different DIC's!
sal      =  33; % Salinity of the sample
tempin   =   7; % Temperature at input conditions
presin   =    100; % Pressure    at input conditions
tempout  =    0; % Temperature at output conditions - doesn't matter in this example
presout  =    0; % Pressure    at output conditions - doesn't matter in this example
sil      =   50; % Concentration of silicate  in the sample (in umol/kg)
po4      =    2; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    10; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

A = CO2SYS_PM(par1,par2,par1type,par2type,sal, ...
    tempin,tempout,presin,presout, ...
    sil,po4,pHscale,k1k2c,kso4c);
PH = A(:,18);
ARAG = A(:,16);
