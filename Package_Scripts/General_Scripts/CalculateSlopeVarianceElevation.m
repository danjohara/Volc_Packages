function [slVarTable,slVarTitles] = CalculateSlopeVarianceElevation(DEM,slp,conts)
% Name: CalculateSlopeVarianceElevation
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
normZ = Zg-minZ;
normZ = normZ./nanmax(normZ(:));

normC = conts-minZ;
normC = normC./nanmax(normC);

slVarTable = [conts(1),normC(1),[1,1,1]*NaN];
slVarTitles = {'True Elevation','Normalized Elevation','Mean Slope','Std Slope','Slope Variance'};
for i = 2:length(normC)
    tSg = Sg;
    tSg(normZ>normC(i)) = NaN;
    tSg(normZ<=normC(i-1)) = NaN;

    meanS = nanmean(tSg(:));
    stdS = nanstd(tSg(:));
    slVar = stdS./meanS;

    slVarTable = [slVarTable;conts(i),normC(i),meanS,stdS,slVar];
end

if slVarTable(end,2) > 1
    slVarTable(end,2) = 1;
end
end