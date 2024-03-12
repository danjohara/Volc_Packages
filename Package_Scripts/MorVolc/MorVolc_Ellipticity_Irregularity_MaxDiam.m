function [ei,ii,contDi,elDi] = MorVolc_Ellipticity_Irregularity_MaxDiam(x,y)
% Name: MorVolc_Ellipticity_Irregularity
% Author: Daniel O'Hara
% Date: 02/24/2021 (mm/dd/yyyy)
% Description: Script to calculate the ellipticity and irregularity index 
%   maximum diameter defined by a set of points (similar to the IDL 
%   version).
%
% Input:
%   x: Array of x-coordinates.
%   y: Array of y-coordinates.
%
% Output:
%   ei: Ellipticity index.
%   ii: Irregularity index.
%   contDi: Boundary dissection index.
%   elDi: Best-fitting ellipse dissection index.

%% Ellipticity Index
[Area,~] = MorVolc_Calculate_Area_Width(x,y);
pd = pdist2([x,y],[x,y]);
L = max(pd(:))/2;
ei = (pi*L^2)/Area;

%% Contour dissection index.
contPerim = 0;
xy = [x,y];
for j = 1:length(x)-1
    contPerim = contPerim + norm(xy(j,:) - xy(j+1,:));
end
contPerim = contPerim + norm(xy(end,:)-xy(1,:));

contDi = contPerim/(2*Area)*sqrt(Area/pi);

%% Best-fitting ellipse dissection index.
A = pi*L^2/ei;
b = A/pi/L;
elPerim = pi*(3*(L+b)-sqrt((3*L+b)*(L+3*b)));

elDi = elPerim/(2*Area)*sqrt(Area/pi);

%% Irregularity Index
ii = contDi-(elDi-1);