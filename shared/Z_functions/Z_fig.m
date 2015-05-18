function [] = Z_fig(fs1)
% 3/2/2011 Parker MacCready
%
% does (hopefully) nice things to figures
%
% you can specify the default fontsize in the argument 'fs1'


if nargin==0; fs1 = 14; end; % fontsize
set(0,'defaultaxesfontsize',fs1);
set(0,'defaulttextfontsize',fs1);
set(0,'defaultaxesfontweight','bold');
set(0,'defaulttextfontweight','bold');
set(0,'defaultaxesfontname','Ariel');

