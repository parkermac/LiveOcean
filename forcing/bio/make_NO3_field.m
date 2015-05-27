function NO3 = make_NO3_field(NO3_method,clm_salt);
% make_NO3_field.m  3/24/2011 Kristen Davis
%
% This code estimates nitrate according to the method chosen by the user
% This is a list of methods and required inputs:
% 1. "PL_Salt" - This is a Piece-wise Linear fit to Salinity and
% requires the input field of salinity
%
% NOTE: Currently only one method!

%-----------------------------------------------------------------------

if NO3_method == 'PL_Salt' % Piece-wise linear fit to salinity (pg 60 in PNWTOX A Notebook)

   % First adjust the NCOM salinity to account for the fresh-bias in order
   % to get the correct amount of nitrate - this correction made by Davis
   % and Siedlecki on 8/18/2011 and updated by KAD on 8/30/12
   salt_corr=zeros(size(clm_salt));
    %salt_corr=clm_salt;
    index=find(clm_salt<=33.45);
    salt_corr(index) = 0.9276*clm_salt(index) + 2.746; % This linear correlation found in comparisons of
    clear index

    index=find(and(clm_salt>33.44,clm_salt<=34.3));
    %salt_corr(index)=clm_salt(index)*0.7828+7.5083; % old salt correction version
    salt_corr(index)=clm_salt(index)*0.6235+12.91;
    clear index

    index=find(clm_salt>34.3);
    salt_corr(index)=clm_salt(index);
    clear index


   % Nitrate-salinity relationship updated by KAD on 8/30/2012
   
    NO3 = zeros(size(salt_corr)); %initializing NO3

    index = find(salt_corr < 31.9);  % Salinity Class 1 (0-31.9 psu)
    NO3(index) = 0;
    clear index

    index = find(and(salt_corr>=31.9, salt_corr<33)); % Salinity Class 2 (31.9-33 psu)
    NO3(index) = 20.15*salt_corr(index) - 642.8;
    clear index

    index = find(and(salt_corr>=33,salt_corr<33.9));  % Salinity Class 3 (33-33.79 psu)
    NO3(index) = 7.7067*salt_corr(index) - 232.1011;
    clear index


    index = find(and(salt_corr>=33.9,salt_corr<34.4));  % Salinity Class 4 (33.79-34.2 psu)
    NO3(index) = 31.6880*salt_corr(index) - 1045.0672;
    clear index


    index = find(and(salt_corr>=34.4,salt_corr<34.43));  % Salinity Class 5 (34.25-34.3 psu)
    NO3(index) = 45;
    clear index

    index = find(salt_corr>=34.43);  % Salinity Class 6 (>34.3)
    NO3(index) = -37.261*salt_corr(index) + 1327.89;
    clear index


   %Old settings
%     index = find(salt_corr >=33.82); % Salinity Class 4 (33.82-34.5 psu)
%     NO3(index) = 34.83*salt_corr(index) - 1148;

   % Now set maximum NO3 to 45 microMolar (found at ~800m depth), based on
   % evidence from historical NO3 data from NODC World Ocean Database
   index = find(NO3 > 45);
   NO3(index) = 45;

   % Also double check that there are no negative values
   clear index
   index = find(NO3 < 0);
   NO3(index) = 0;

else disp(['NO3 estimation method not recognized.'])
end