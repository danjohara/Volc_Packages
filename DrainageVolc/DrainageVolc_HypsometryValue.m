function [areaHyps,valHyps] = DrainageVolc_HypsometryValue(valueGrid,iter,normalize)
% Name: HypsometryValue
% Author: Daniel O'Hara
% Date: 02/04/2021 (mm/dd/yyyy)
% Description: Script to calculate the hypsometry (CDF) of a given
%   distribution of values, either normalized or non-normalized.
%
% Input:
%   valueGrid: Grid of values to analyze.
%   iter: Iteration value for hypsometry calculation.
%   normalize: Flag to normalize results.
%
% Output:
%   areaHyps: Normalized area of values.
%   valHyps: Distribution of values.

totCells = sum(sum(~isnan(valueGrid)));
offset = min(valueGrid(:));
valueGrid = valueGrid - offset;
coeff = max(valueGrid(:));
valueGrid = valueGrid./coeff;

valHyps = 0:iter:1;
areaHyps = zeros(size(valHyps));
for i = 1:length(valHyps)
    areaHyps(i) = sum(sum(valueGrid<=valHyps(i)))./totCells;
end
if ~normalize
    valHyps = valHyps*coeff + offset;
end