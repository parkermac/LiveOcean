function RH = Z_wmo_RH(P,T,Q)
% 5/21/2011 Nick Lederer, modified by Parker MacCready
%
% this converts mixing ratio (kg kg-1) which is the usual mm5 output, into
% relative humidity (%) which is what ROMS expects
%
% INPUT:
% P in hectaPascal or millibar
% T in Celcius
% Q in kg kg-1
%
% OUTPUT:
% RH in percent

% all equations come from Chapter 4 of
% http://www.wmo.int/pages/prog/www/IMOP/publications/CIMO-Guide/
% CIMO_Guide-7th_Edition-2008.html
%
% WMO GUIDE TO METEOROLOGICAL INSTRUMENTS AND METHODS OF OBSERVATION
%          WMO-No. 8 (Seventh edition) (6 August 2008)

e_prime = Q.*P./(0.62198+Q); % WHO equation 4.A.6
fp = 1.0016 + 3.15e-6 * P - 0.074./P; % from Annex 4.B
ew = 6.112*exp(17.62.*T./(243.12 + T)); % from Annex 4.B
ew_prime = fp.*ew; % from Annex 4.B
RH = 100 * e_prime./ew_prime; % from Annex 4.B
