function oxygen = make_oxy_field(NO3_method,salt)
%
% make_oxy_field.m
%
% This code estimates oxygen according to the method chosen by the user.
% 
% Supported methods:
%
%  PL_salt = piecewise-linear regression on salt
%
% Recoded 9/16/2015 by PM to clean up and use of new regressions
% from Ryan McCabe of 8/2015.

%-----------------------------------------------------------------------

if strcmp(NO3_method,'PL_Salt') && (nanmax(salt(:)) < 35)
    
    % Salinity vs. oxygen [uM], Ryan McCabe 8/2015
    % oxygen = mm*salt + bb;
        
    mm = zeros(size(salt));
    bb = zeros(size(salt));
    
    ind = (salt < 32.167);
    mm(ind) = 0;
    bb(ind) = 300;
    ind = ((salt >= 32.167) & (salt < 33.849));
    mm(ind) = -113.9481;
    bb(ind) = 3965.3897;
    ind = ((salt >= 33.849) & (salt < 34.131));
    mm(ind) = -278.3006;
    bb(ind) = 9528.5742;
    ind = ((salt >= 34.131) & (salt < 34.29));
    mm(ind) = -127.2707;
    bb(ind) = 4373.7895;
    ind = ((salt >= 34.29) & (salt < 34.478));
    mm(ind) = 34.7556;
    bb(ind) = -1182.0779;
    ind = ((salt >= 34.478) & (salt < 35));
    mm(ind) = 401.7916;
    bb(ind) = -13836.8132;
        
    oxygen = mm.*salt + bb;
    
    % Limit values.
    oxygen(oxygen > 450) = 450;
    oxygen(oxygen < 0) = 0;
    
else
    
    disp('Estimation method not recognized, or salt out of range.')
    
end