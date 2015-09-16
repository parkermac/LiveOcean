function [DIC, TAlk] = make_DIC_field(NO3_method,salt)
%
% make_DIC_field.m
%
% This code estimates DIC & TAlk according to the method chosen by the user.
% 
% Supported methods:
%
%  PL_salt = piecewise-linear regression on salt
%
% Recoded 9/16/2015 by PM to clean up and use of new regressions
% from Ryan McCabe of 8/2015.

%-----------------------------------------------------------------------

if strcmp(NO3_method,'PL_Salt') && (nanmax(salt(:)) < 35)
    
    % Salinity vs. TIC [uM]
    % TIC = mm*salt + bb;
        
    mm = zeros(size(salt));
    bb = zeros(size(salt));
    
    ind = (salt < 31.887);
    mm(ind) = 27.7967;
    bb(ind) = 1112.2027;
    ind = ((salt >= 31.887) & (salt < 33.926));
    mm(ind) = 147.002;
    bb(ind) = -2688.8534;
    ind = ((salt >= 33.926) & (salt < 34.197));
    mm(ind) = 352.9123;
    bb(ind) = -9674.5448;
    ind = ((salt >= 34.197) & (salt < 34.504));
    mm(ind) = 195.638;
    bb(ind) = -4296.2223;
    ind = ((salt >= 34.504) & (salt < 35));
    mm(ind) = -12.7457;
    bb(ind) = 2893.77;
            
    DIC = mm.*salt + bb;

    % Salinity vs. TAlk [uM]
    % TAlk = mm*salt + bb;
    
    ind = (salt < 31.477);
    mm(ind) = 37.0543;
    bb(ind) = 1031.0726;
    ind = ((salt >= 31.477) & (salt < 33.915));
    mm(ind) = 48.5821;
    bb(ind) = 668.2143;
    ind = ((salt >= 33.915) & (salt < 35));
    mm(ind) = 246.2214;
    bb(ind) = -6034.6841;
            
    TAlk = mm.*salt + bb;
        
else
    
    disp('Estimation method not recognized, or salt out of range.')
    
end