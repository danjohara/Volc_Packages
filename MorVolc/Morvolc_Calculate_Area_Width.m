function [Area,Width] = Morvolc_Calculate_Area_Width(x,y)
% Name: Morvolc_Calculate_Area_Width
% Author: Daniel O'Hara
% Date: 02/24/2021 (mm/dd/yyyy)
% Description: Script to calculate the area and equivalant circle diameter
%   of a given boundary.
%
% Input:
%   x: X-coordinates of boundary.
%   y: Y-coordinates of boundary.
%
% Output:
%   Area: Area of the region encompassed by boundary.
%   Width: Diameter of circle that has equivalant area.

Area = polyarea(x,y);
Width = 2*sqrt(Area/pi);
end