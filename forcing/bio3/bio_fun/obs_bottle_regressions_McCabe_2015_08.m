% obs_bottle_regressions.m
%
% Simple storage for piecewise linear regression coefficients that were
% constructed from PMEL bottle data.
%
% R. McCabe 08/2015


%% Salinity vs. TA [uM]
% TA = mm*salt + bb;
if salt < 31.477
    mm = 37.0543;
    bb = 1031.0726;
elseif ((salt >= 31.477) && (salt < 33.915))
    mm = 48.5821;
    bb = 668.2143;
elseif ((salt >= 33.915) && (salt < 35))
    mm = 246.2214;
    bb = -6034.6841;
end


%% Salinity vs. TIC [uM]
% TIC = mm*salt + bb;
if salt < 31.887
    mm = 27.7967;
    bb = 1112.2027;
elseif ((salt >= 31.887) && (salt < 33.926))
    mm = 147.002;
    bb = -2688.8534;
elseif ((salt >= 33.926) && (salt < 34.197))
    mm = 352.9123;
    bb = -9674.5448;
elseif ((salt >= 34.197) && (salt < 34.504))
    mm = 195.638;
    bb = -4296.2223;
elseif ((salt >= 34.504) && (salt < 35))
    mm = -12.7457;
    bb = 2893.77;
end


%% Salinity vs. Oxygen [uM]
% Oxygen = mm*salt + bb;
if salt < 32.167
    mm = 0;
    bb = 300;
elseif ((salt >= 32.167) && (salt < 33.849))
    mm = -113.9481;
    bb = 3965.3897;
elseif ((salt >= 33.849) && (salt < 34.131))
    mm = -278.3006;
    bb = 9528.5742;
elseif ((salt >= 34.131) && (salt < 34.29))
    mm = -127.2707;
    bb = 4373.7895;
elseif ((salt >= 34.29) && (salt < 34.478))
    mm = 34.7556;
    bb = -1182.0779;
elseif ((salt >= 34.478) && (salt < 35))
    mm = 401.7916;
    bb = -13836.8132;
end


%% Salinity vs. NO3 [uM]
% NO3 = mm*salt + bb;
if salt < 31.898
    mm = 0;
    bb = 0;
elseif ((salt >= 31.898) && (salt < 33.791))
    mm = 16.3958;
    bb = -522.989;
elseif ((salt >= 33.791) && (salt < 34.202))
    mm = 29.6973;
    bb = -972.4545;
elseif ((salt >= 34.202) && (salt < 34.482))
    mm = 8.0773;
    bb = -233.0007;
elseif ((salt >= 34.482) && (salt < 35))
    mm = -28.6251;
    bb = 1032.5686;
end


