function worker(tempdir)

% worker to call CO2SYS.m

a = load([tempdir,'worker_input.mat']);

A = CO2SYS(a.alkalinity(:), a.TIC(:), 1, 2, a.salt(:), a.temp(:), a.temp(:), a.pres(:), a.pres(:), 50, 2, 1, 10, 1);

ph = A(:,18);
om = A(:,31);

PH = reshape(ph, size(a.salt));
OM = reshape(om, size(a.salt));

%PH = real(PH);
%OM = real(OM);

save([tempdir,'PH.mat'], 'PH');
save([tempdir,'OM.mat'], 'OM');