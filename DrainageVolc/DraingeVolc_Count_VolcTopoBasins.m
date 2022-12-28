function [normTopo,numDB] = DraingeVolc_Count_VolcTopoBasins(Z,DB,iter)
% Name: CountVolcTopoBasins
% Author: Daniel O'Hara
% Date: 02/04/2021 (mm/dd/yyyy)
% Description: Script to count the number of basins that exist with 
%   increasing elevation of a DEM.
%
% Input:
%   Z: Grid of elevations.
%   DB: Grid of drainage basins.
%   iter: Iteration value for increasing normalized topography (equivalant 
%       to hypsometry calculation).
%
% Output:
%   normTopo: Normalized elevation values.
%   numDB: Number of basins, corresponding to normalized elevation.

Z = Z-min(Z(:));
Z = Z./max(Z(:));

normTopo = 0:iter:1;
numDB = zeros(size(normTopo));
for i = 1:length(normTopo)
    t = Z >= normTopo(i);
    
    dd = DB(t==1);
    numDB(i) = numel(unique(dd));
end