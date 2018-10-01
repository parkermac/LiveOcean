% code to test CO2SYS

% Assume in mean in situ, and out means at the surface (changes pressure)
PAR1 = 2400; %  (some unit) : scalar or vector of size n
PAR2 = 2200; %  (some unit) : scalar or vector of size n
PAR1TYPE = 1; %      () : scalar or vector of size n (*)
PAR2TYPE = 2; %     () : scalar or vector of size n (*)
SAL = 35; %            () : scalar or vector of size n
TEMPIN = 8; %  (degr. C) : scalar or vector of size n
TEMPOUT = 8; % (degr. C) : scalar or vector of size n
PRESIN = 4000; %     (dbar) : scalar or vector of size n
PRESOUT = 0; %   (dbar) : scalar or vector of size n
SI = 50; %    (umol/kgSW) : scalar or vector of size n
PO4 = 2; %   (umol/kgSW) : scalar or vector of size n
pHSCALEIN = 1; %        : scalar or vector of size n (**)
K1K2CONSTANTS = 10; %     : scalar or vector of size n (***)
KSO4CONSTANTS = 1; %    : scalar or vector of size n (****)

[RESULT,HEADERS,NICEHEADERS]=CO2SYS_PM(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
    SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,...
    K1K2CONSTANTS,KSO4CONSTANTS);

np = size(RESULT,2);
for ii = 1:np
    fprintf('%0.2f =  %s\n',RESULT(ii),NICEHEADERS{ii})
end

