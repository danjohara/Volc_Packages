function [slVarTable,slVarTitles] = DrainageVolc_CalculateSlopeVarianceElevation(DEM,slp,contIter)
% Name: DrainageVolc_CalculateSlopeVarianceElevation
% Author: Daniel O'Hara
% Date: 08/18/2023 (mm/dd/yyyy)
% Description: Script to calculate the slope variance (std / mean) of
% elevation bins.
%
% Input:
%   DEM: GRIDobj of elevations.
%   slp: GRIDobj of slope.
%   slpSize: Contour interval for bins. Negative values between -1 and 0
%       are treated as a percentage of the total edifice height.
%
% Output:
%   slVarTable: Array of values. Columns are elevation bin, normalized
%       elevation bin, mean slope, std slope, and slope variance.
%   slVarTitles: Cell array describing slVarTable columns.

[Zg,~,~] = GRIDobj2mat(DEM);
[Sg,~,~] = GRIDobj2mat(slp);

minZ = nanmin(Zg(:));
Zg = Zg - minZ;
if contIter < 0 && contIter >= -1
    contIter = range(Zg(:))*abs(contIter);
end

slVarTable = [minZ,0,[1,1,1]*NaN];
slVarTitles = {'True Elevation','Normalized Elevation','Mean Slope','Std Slope','Slope Variance'};
for i = contIter:contIter:max(Zg(:))+contIter
    tSg = Sg;
    tSg(Zg > i) = NaN;
    tSg(Zg < i-contIter) = NaN;

    meanS = nanmean(tSg(:));
    stdS = nanstd(tSg(:));
    slVar = stdS./meanS;

    slVarTable = [slVarTable;i+minZ,i/range(Zg(:)),meanS,stdS,slVar];
end

if slVarTable(end,2) > 1
    slVarTable = slVarTable(1:end-1,:);
end