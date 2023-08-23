function [normTopo,numDB,numDB2] = DraingeVolc_Count_VolcTopoBasins(Z,DB,A,iter,areaThreshold)
% Name: CountVolcTopoBasins
% Author: Daniel O'Hara
% Date: 02/04/2021 (mm/dd/yyyy)
% Description: Script to count the number of basins that exist with 
%   increasing elevation of a DEM.
%
% Input:
%   Z: Grid of elevations.
%   DB: Grid of drainage basins.
%   A: GRIDobj of maximum basin drainage areas.
%   iter: Iteration value for increasing normalized topography (equivalant 
%       to hypsometry calculation).
%   areaThreshold: Drainage area threshold for secondary basin count.
%
% Output:
%   normTopo: Normalized elevation values.
%   numDB: Number of basins, corresponding to normalized elevation.

Z = Z-min(Z(:));
Z = Z./max(Z(:));

DB2 = DB;
DB2(A<areaThreshold) = NaN;

normTopo = 0:iter:1;
numDB = zeros(size(normTopo));
numDB2 = numDB;
for i = 1:length(normTopo)
    t = Z >= normTopo(i);
    
    dd = DB(t==1);
    tmp = unique(dd);
    tmp(isnan(tmp)) = [];
    numDB(i) = numel(tmp);

    dd = DB2(t==1);
    tmp = unique(dd);
    tmp(isnan(tmp)) = [];
    numDB2(i) = numel(tmp);
end