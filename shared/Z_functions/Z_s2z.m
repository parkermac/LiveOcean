function [z_rho,z_w] = Z_s2z(h,zeta,S)
% useage: [z_rho,z_w] = Z_s2z(h,zeta,S)
% Z_s2z.m  12/19/2012  Parker MacCready & Sarah Giddings
%
% This gives up to 3D arrays of z position (m, positive up)
% packed from the bottom to the top, for a given bottom depth (h)
% and surface height
% Note that the inputs "h" and "zeta" are matrices [LxM],
% and "S" is a structure
% created by "Z_get_basic_info.m" which contains all the s-coordinate info.
% z_rho is at box mid points, and z_w is at box top or bottom edges
%
% IMPORTANT:  SNG Mar 2011 extended this code to include Vtranform = 2 (see 
% https://www.myroms.org/wiki/index.php/Vertical_S-coordinate for transform
% information. Also vectorized for speed.
% Also note that if you choose either Vtransform = 1 or 2 and set hc = 0
% then you can use roms_z.m which is slightly faster.

%put all parameters in same matrix format
[M,N] = size(h);
[L] = length(S.s_rho);
s_rho = repmat(S.s_rho,[1 M N]);
Cs_r = repmat(S.Cs_r,[1 M N]);
s_w = repmat(S.s_w,[1 M N]);
Cs_w = repmat(S.Cs_w,[1 M N]);
Hc_r = repmat(S.hc,[L M N]);
H_r = repmat(reshape(h,[1 M N]),[L 1 1]);
Zeta_r = repmat(reshape(zeta,[1 M N]),[L 1 1]);
Hc_w = repmat(S.hc,[L+1 M N]);
H_w = repmat(reshape(h,[1 M N]),[L+1 1 1]);
Zeta_w = repmat(reshape(zeta,[1 M N]),[L+1 1 1]);

%if hc = 0, eqns are simpler (hence why roms_z is faster) and are the same
%regardless of if Vtransform = 1 or 2
if S.hc == 0
    z_rho = H_r.*Cs_r+Zeta_r+Zeta_r.*(Cs_r);
    z_w = H_w.*Cs_w+Zeta_w+Zeta_w.*(Cs_w);
else
    %eqn. for Vtransform = 1
    if S.Vtransform == 1;
        zr0 = (s_rho-Cs_r).*Hc_r+Cs_r.*H_r;
        z_rho = zr0+Zeta_r.*(1+(zr0./H_r));
        zw0 = (s_w-Cs_w).*Hc_w+Cs_w.*H_w;
        z_w = zw0+Zeta_w.*(1+(zw0./H_w));
    %eqn. for Vtransform = 2
    elseif S.Vtransform == 2;
        zr0 = ((s_rho.*Hc_r)+(Cs_r.*H_r))./(Hc_r+H_r);
        z_rho = Zeta_r + (Zeta_r+H_r).*zr0;
        zw0 = ((s_w.*Hc_w)+(Cs_w.*H_w))./(Hc_w+H_w);
        z_w = Zeta_w + (Zeta_w+H_w).*zw0;
    end
end
