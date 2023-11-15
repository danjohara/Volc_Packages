function totalDD = DrainageVolc_Collect_Total_DrainageDensity(Z,S,A,DB,dxy)
% Name: Collect_Total_DrainageDensity
% Author: Daniel O'Hara
% Date: 04/2/2021 (mm/dd/yyyy)
% Description: Script to calculate the cumulative drainage density over the
%   entire landform.
%
% Input:
%   Z: Matrix of elevations.
%   S: STREAMobj of all edifice streams.
%   A: GRIDobj of accumulating drainage areas.
%   DB: GRIDboj of drainage basins.
%   dxy: Grid spacing.
%
% Output:
%   totalDD: Drainage density, using the cumulative length of channels and
%       the entire drainage area of the edifice.

%% Setup
[DBg,x,y] = GRIDobj2mat(DB);
DBg(isnan(Z)) = NaN;
[X,Y] = meshgrid(x,y);
dbi = unique(DBg(:));
dbi(isnan(dbi)) = [];

[Ag,~,~] = GRIDobj2mat(A);
Ag(isnan(Z)) = NaN;
%% Collect Distances and Make Grid
streamDistances = distance(S,'accumdownstream');
sD = STREAMobj2GRIDobj(S,streamDistances);
[sDg,~,~] = GRIDobj2mat(sD);

%% Get Total Stream Lengths and All Areas
totalS = 0;
totalS2 = 0;
totalA = 0;
for i = 1:length(dbi)
    tmpAg = Ag;
    tmpAg(DBg~=dbi(i)) = NaN;
    largestA = nanmax(tmpAg(:));
    if ~isnan(largestA)
        totalA = totalA + largestA;
    end
    
    tmpSDg = sDg;
    tmpSDg(DBg~=dbi(i)) = NaN;
    largestL = max(tmpSDg(:));
    if ~isnan(largestL)
        totalS = totalS+largestL;
    end
end

%% Calculate Density
totalDD = totalS/(totalA*dxy^2);
end