function [] = Z_info(basename,tt,T,poschar)
% 3/2/2011 Parker MacCready
%
% this adds informational text from a ROMS history file to an existing plot
%
% poschar is the position 'ul','ur','ll','lr'

[xt,yt]=Z_lab(poschar);
aa = axis; Dlat = aa(4)-aa(3);
dy = Dlat/30;

% this sets the horizontalalignment character
hal = poschar(2);

% this sets the vertical offset
switch poschar(1)
    case 'u'
        voff = -3;
    case 'l'
        voff = 0;
end

% put the text on the figure
text(xt,yt + (voff+3)*dy,[datestr(T.time_datenum,1),' '], ...
    'horizontalalignment',hal);
text(xt,yt + (voff+2)*dy,[datestr(T.time_datenum,15),' GMT '], ...
    'horizontalalignment',hal);
text(xt,yt + (voff+1)*dy,[strrep(basename,'_',' '),' '], ...
    'horizontalalignment',hal);
text(xt,yt + (voff+0)*dy,['save ',num2str(tt),' '], ...
    'horizontalalignment',hal);
