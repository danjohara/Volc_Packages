function [BCount_Table,BCount_Titles,BIDs,BIDs_AboveA,R_DEM,RNorm_DEM] = DrainageVolc_CalculateRadialDistanceNumbers(DEM,DB,A,distIter,areaThreshold)
% Name: DrainageVolc_CalculateRadialDistanceNumbers
% Author: Daniel O'Hara
% Date: 06/30/2021 (mm/dd/yyyy)
% Description: Script to calculate the number of basins (and number of 
%    basins above a channelization threshold) at intervaled distances from 
%    the edifice's peak. Distance values are normalized relative to the
%    maxiimum edifice radius.
%
% Input:
%   DEM: GRIDobj of elevations.
%   DB: GRIDobj of drainage basins.
%   A: GRIDobj of maximum basin drainage areas.
%   distIter: Interval for distance counts. If the value is between -1 and 
%       0, distance is treated as a percentage of the maximum radius.                                                                             
%   areaThreshold: Drainage area threshold for secondary basin count.
%
% Output:
%   BCount_Table: Table of results. Columns are bin distance, normalized
%       bin distance, all basin count, and basin count of only those above
%       areaThreshold.
%   BCount_Titles: Cell array of table columns.
%   BIDs: Basin ID's at each interval.
%   BIDs_AboveA: Basin ID's above the channelization threshold at each
%       interval.
%   R_DEM: GRIDobj of radial distances from the edifice's peak.
%   RNorm_DEM: GRIDobj of raidal distances from the edifice's peak,
%       normalized by the maximum radius.

%% Get radial distances
[Z,x,y] = GRIDobj2mat(DEM);
[X,Y] = meshgrid(x,y);

[i,j] = find(Z==max(Z(:)),1);
R = sqrt((X-X(i,j)).^2 + (Y-Y(i,j)).^2);
R(isnan(Z)) = NaN;

%% Normalize distances
R_Norm = R./max(R(:));

%% Count Basins
[DBg,~,~] = GRIDobj2mat(DB);
[Ag,~,~] = GRIDobj2mat(A);

% Perpare secondary count
DBg2 = DBg;
if ~isnan(areaThreshold)
    DBg2(Ag<areaThreshold) = NaN;
end

% Create normalize values
if distIter < 0 && distIter >= -1
    normDistIter = abs(distIter);
    distIter = normDistIter*max(R(:));
elseif distIter > 0
    normDistIter = distIter./max(R(:));
end

% Make bins
normDistBins = normDistIter:normDistIter:1+normDistIter;
normDistBins(normDistBins>1) = 1;
normDistBins = unique(normDistBins);
distBins = normDistBins*max(R(:));

BCounts = zeros(size(normDistBins));
BCounts_AboveA = BCounts;
BIDs = cell(size(normDistBins));
BIDs_AboveA = BIDs;

for i = 1:length(normDistBins)
    %% All basins
    tmp = DBg;
    tmp(R_Norm>normDistBins(i)) = NaN;
    uniDB = unique(tmp(:));
    uniDB(isnan(uniDB)) = [];

    BCounts(i) = length(uniDB);
    BIDs{i} = uniDB;

    %% Basins above the threshold
    tmp = DBg2;
    tmp(R_Norm>normDistBins(i)) = NaN;
    uniDB = unique(tmp(:));
    uniDB(isnan(uniDB)) = [];

    BCounts_AboveA(i) = length(uniDB);
    BIDs_AboveA{i} = uniDB;
end

%% Make Table
BCount_Table = [distBins',normDistBins',BCounts',BCounts_AboveA'];
BCount_Titles = {'Distance','Normalized Distance','Number of Basins','Number of Basins Above Threshold'};
R_DEM = GRIDobj(X,Y,R);
RNorm_DEM = GRIDobj(X,Y,R_Norm);
end
