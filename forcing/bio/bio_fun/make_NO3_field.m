function NO3 = make_NO3_field(NO3_method,salt);
%
% make_NO3_field.m
%
% This code estimates nitrate according to the method chosen by the user.
% 
% Supported methods:
%
%  PL_salt = piecewise-linear regression on salt
%
% Recoded 9/16/2015 by PM to clean up and use of new regressions
% from Ryan McCabe of 8/2015.

%-----------------------------------------------------------------------

if strcmp(NO3_method,'PL_Salt') && (nanmax(salt(:)) < 35)
    
    % Salinity vs. NO3 [uM], Ryan McCabe 8/2015
    % NO3 = mm*salt + bb;
        
    mm = zeros(size(salt));
    bb = zeros(size(salt));
    
    ind = (salt < 31.898);
    mm(ind) = 0;
    bb(ind) = 0;
    ind = ((salt >= 31.898) & (salt < 33.791));
    mm(ind) = 16.3958;
    bb(ind) = -522.989;
    ind = ((salt >= 33.791) & (salt < 34.202));
    mm(ind) = 29.6973;
    bb(ind) = -972.4545;
    ind = ((salt >= 34.202) & (salt < 34.482));
    mm(ind) = 8.0773;
    bb(ind) = -233.0007;
    ind = ((salt >= 34.482) & (salt < 35));
    mm(ind) = -28.6251;
    bb(ind) = 1032.5686;
    
    NO3 = mm.*salt + bb;
    
    % Set maximum NO3 to 45 microMolar (found at ~800m depth), based on
    % evidence from historical NO3 data in NODC World Ocean Database.
    NO3(NO3 > 45) = 45;
    
    % Ensure that there are no negative values.
    NO3(NO3 < 0) = 0;
    
else
    
    disp('Estimation method not recognized, or salt out of range.')
    
end