% test reading lo_info

addpath('/Users/pm7/Documents/LiveOcean/alpha')

% find the path to alpha (assumes math has been added
% by the calling function
alp0 = which('Lstart');
alp = alp0(1:end-8);

fid = fopen([alp,'lo_info.csv']);
while ~feof(fid)
    line = fgetl(fid);
    ic = find(line == ',');
    name = line(1:ic-1);
    value = line(ic+1:end);
    Ldir.(name) = value;
end
fclose(fid);