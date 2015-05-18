function [xt,yt] = Z_lab(poschar)
% 3/2/2011 Parker MacCready
%
% Returns the position for a text label
%
% Say you have a plot up, and "aa = axis"
% then use poschar = 'll', 'ul', 'ur', or 'lr'
% (meaning "lower left' and etc.)
% to define the position you want

aa = axis;

if ischar(poschar) & length(poschar)==2
   switch poschar
   case 'll'
      xt = aa(1) + 0.05 * (aa(2)-aa(1));
      yt = aa(3) + 0.10 * (aa(4)-aa(3));
   case 'lr'
      xt = aa(2) - 0.05 * (aa(2)-aa(1));
      yt = aa(3) + 0.10 * (aa(4)-aa(3));
   case 'ul'
      xt = aa(1) + 0.05 * (aa(2)-aa(1));
      yt = aa(4) - 0.10 * (aa(4)-aa(3));
   case 'ur'
      xt = aa(2) - 0.05 * (aa(2)-aa(1));
      yt = aa(4) - 0.10 * (aa(4)-aa(3));
   end
else
   disp('Need to input a a different poschar')
   return
end

