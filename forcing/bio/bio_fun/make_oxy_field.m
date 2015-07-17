function oxygen = make_oxy_field(NO3_method,clm_salt)
% make_oxy_field.m  4/1/2011 Samantha Siedlecki after Kristen Davis'
% Nitrate function
%
% This code estimates nitrate according to the method chosen by the user
% This is a list of methods and required inputs:
% 1. "PL_Salt" - This is a Piece-wise Linear fit to Salinity and
% requires the input field of salinity
%
% NOTE: Currently only one method!

if strcmp(NO3_method,'PL_Salt') % Piece-wise linear fit to salinity (pg 60 in PNWTOX A Notebook)
    salt_corr=clm_salt; % legacy code from when we needed to correct NCOM salt bias
    
    oxygen = zeros(size(salt_corr)); %initializing oxygen
    
    oxygen= (salt_corr)*-129.23 + 4482.4;
    
    index=find(and(salt_corr>=33.9,salt_corr<34.5));
    oxygen(index)=511.42*salt_corr(index).^2-35084*salt_corr(index)+601730;
    clear index
    
    index = find(salt_corr>=34.5);
    oxygen(index) = 4e-100*exp(6.7368*(salt_corr(index)));
    clear index
    
    % Now set minimum of 0 for oxygen s
    index = find(oxygen < 0);
    oxygen(index) = 0;
    % Now set maximum of 450 for oxygen
    index = find(oxygen > 450);
    oxygen(index) = 450;
    
else disp(['oxygen estimation method not recognized.'])
end