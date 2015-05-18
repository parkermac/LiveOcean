function [z_w] = Z_s2z_w(h,zeta,S)
% Z_s2z_w.m  1/30/2013  Parker MacCready & Sarah Giddings
%
% USEAGE:
% [z_w] = Z_s2z_w(h,zeta,S)
%
% Just like Z_s2z.m but only returns z_w

%put all parameters in same matrix format
[M,N] = size(h);
[L] = length(S.s_rho);
s_w = repmat(S.s_w,[1 M N]);
Cs_w = repmat(S.Cs_w,[1 M N]);
Hc_w = repmat(S.hc,[L+1 M N]);
H_w = repmat(reshape(h,[1 M N]),[L+1 1 1]);
Zeta_w = repmat(reshape(zeta,[1 M N]),[L+1 1 1]);

%if hc = 0, eqns are simpler (hence why roms_z is faster) and are the same
%regardless of if Vtransform = 1 or 2
if S.hc == 0
    z_w = H_w.*Cs_w+Zeta_w+Zeta_w.*(Cs_w);
else
    %eqn. for Vtransform = 1
    if S.Vtransform == 1;
        zw0 = (s_w-Cs_w).*Hc_w+Cs_w.*H_w;
        z_w = zw0+Zeta_w.*(1+(zw0./H_w));
    %eqn. for Vtransform = 2
    elseif S.Vtransform == 2;
        zw0 = ((s_w.*Hc_w)+(Cs_w.*H_w))./(Hc_w+H_w);
        z_w = Zeta_w + (Zeta_w+H_w).*zw0;
    end
end
