function [] = Z_dar
% 3/2/2011 Parker MacCready
%
% fixes aspect ratio of Cartesian plots so that the aspect ratio is correct
% based on the average latitude of the plot

aa = axis;
meanlat = mean(aa(3:4));
dar = [1/cos(pi*meanlat/180) 1 1]; % Cartesian scaling
set(gca,'dataaspectratio',dar);
