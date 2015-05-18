function S = Z_scoord(theta_s,theta_b,tcline,hmin,N,Vtransform,Vstretching)
% Z_scoord.m  5/21/2007  Parker MacCready
%
% this creates the structure S, which would be used for example by
% Z_s2z.m, given basic grid parameters
%
% IMPORTANT: note that this is the Song & Haidvogel 1994 stretching
% function, ROMS Vstretching = 1. If one chooses a different stretching
% function (as of March 2011 Vstretching can be 1, 2, or 3) then this code
% must be updated to include the proper stretching transformation!
%
% edited by DAS to include more things in S stucture
% edited by SNG March 2011 to include all of the current available ROMS
% stretching functions, 1-4 see: 
% https://www.myroms.org/wiki/index.php/Vertical_S-coordinate#Vertical_Stretching_Functions
%-----------------------------------------------------------------------

if Vtransform == 1
    hc = min(hmin,tcline);
elseif Vtransform == 2
    hc = tcline;
end

sc_r = ([-(N-1):0]-0.5)/N;
sc_w = [-N:0]/N;

if Vstretching == 1
    if theta_s ~= 0
        cff1 = 1/sinh(theta_s);
        cff2 = 0.5/tanh(0.5*theta_s);
        Cs_r = (1-theta_b)*cff1*sinh(theta_s*sc_r) ...
            + theta_b*( cff2*tanh(theta_s*(sc_r + 0.5)) - 0.5 );
        Cs_w = (1-theta_b)*cff1*sinh(theta_s*sc_w) ...
            + theta_b*( cff2*tanh(theta_s*(sc_w + 0.5)) - 0.5 );
    else
        Cs_r = sc_r;
        Cs_w = sc_w;
    end
elseif Vstretching == 2
    alpha = 1; beta = 1;
    if theta_s~=0 && theta_b~=0
        Csur = (1-cosh(theta_s.*sc_r))./(cosh(theta_s)-1);
        Cbot = ((sinh(theta_b.*(sc_r+1)))./(sinh(theta_b)))-1;
        u = ((sc_r+1).^alpha).*(1+(alpha/beta)*(1-((sc_r+1).^beta)));
        Cs_r = u.*Csur+(1-u).*Cbot;
        Csur_w = (1-cosh(theta_s.*sc_w))./(cosh(theta_s)-1);
        Cbot_w = ((sinh(theta_b.*(sc_w+1)))./(sinh(theta_b)))-1;
        u_w = ((sc_w+1).^alpha).*(1+(alpha/beta)*(1-((sc_w+1).^beta)));
        Cs_w = u_w.*Csur_w+(1-u_w).*Cbot_w;
    else
        Cs_r = sc_r;
        Cs_w = sc_w;
    end
elseif Vstretching == 3
    %Geyer function for high bbl resolution in shallow applications
    gamma = 3;
    Csur = -(log(cosh(gamma.*abs(sc_r).^theta_s)))./log(cosh(gamma));
    Cbot = ((log(cosh(gamma.*(sc_r+1).^theta_b)))./log(cosh(gamma)))-1;
    mu = 0.5*(1-tanh(gamma*(sc_r+0.5)));
    Cs_r = mu.*Cbot+(1-mu).*Csur;
    Csur_w = -(log(cosh(gamma.*abs(sc_w).^theta_s)))./log(cosh(gamma));
    Cbot_w = ((log(cosh(gamma.*(sc_w+1).^theta_b)))./log(cosh(gamma)))-1;
    mu_w = 0.5*(1-tanh(gamma*(sc_w+0.5)));
    Cs_w = mu_w.*Cbot_w+(1-mu_w).*Csur_w;    
elseif Vstretching == 4
    %newest ROMS default as of March 2011 (theta_s between 0 and 10,
    % theta_b between 0 and 4)
    if theta_s>0
        Cs_r = (1-cosh(theta_s.*sc_r))./(cosh(theta_s)-1);
        Cs_w = (1-cosh(theta_s.*sc_w))./(cosh(theta_s)-1);
    elseif theta_s<=0
        Cs_r = -(sc_r.^2);
        Cs_w = -(sc_w.^2);
    end
    if theta_b > 0
        Cs_r = (exp(theta_b.*Cs_r)-1)./(1-exp(-theta_b));
        Cs_w = (exp(theta_b.*Cs_w)-1)./(1-exp(-theta_b));
    end
end


S.s_rho = sc_r';
S.s_w = sc_w';
S.hc = hc;
S.Cs_r = Cs_r';
S.Cs_w = Cs_w';
S.N = N;
S.theta_b = theta_b;
S.hmin = hmin;
S.theta_s = theta_s;
S.tcline = tcline; 
S.Vtransform = Vtransform;
S.Vstretching = Vstretching;