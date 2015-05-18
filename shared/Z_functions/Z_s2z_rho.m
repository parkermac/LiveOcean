function [z_rho] = Z_s2z_rho(h,zeta,S)
% Z_s2z_rho.m  1/30/2013  Parker MacCready & Sarah Giddings
%
% Just like Z_s2z.m but only returns z_rho

%put all parameters in same matrix format
[M,N] = size(h);
[L] = length(S.s_rho);
s_rho = repmat(S.s_rho,[1 M N]);
Cs_r = repmat(S.Cs_r,[1 M N]);
Hc_r = repmat(S.hc,[L M N]);
H_r = repmat(reshape(h,[1 M N]),[L 1 1]);
Zeta_r = repmat(reshape(zeta,[1 M N]),[L 1 1]);

%if hc = 0, eqns are simpler (hence why roms_z is faster) and are the same
%regardless of if Vtransform = 1 or 2
if S.hc == 0
    z_rho = H_r.*Cs_r+Zeta_r+Zeta_r.*(Cs_r);
else
    %eqn. for Vtransform = 1
    if S.Vtransform == 1;
        zr0 = (s_rho-Cs_r).*Hc_r+Cs_r.*H_r;
        z_rho = zr0+Zeta_r.*(1+(zr0./H_r));
    %eqn. for Vtransform = 2
    elseif S.Vtransform == 2;
        zr0 = ((s_rho.*Hc_r)+(Cs_r.*H_r))./(Hc_r+H_r);
        z_rho = Zeta_r + (Zeta_r+H_r).*zr0;
    end
end
