function [Ss,S_Stats] = DrainageVolc_Collect_Stream_Concavity(bIAL,DEM,A,DB,FD,basinIDs,channelThreshold)
% Name: SeparateKLargestStreams
% Author: Daniel O'Hara
% Date: 05/11/2021 (mm/dd/yyyy)
% Description: Script to calculate the concavity indexes of the largest 
%   basins. 
%
% Input:
%   bIAL: Array of basin statistics from Collect_DrainageBasin_Stats_wCross.
%   DEM: GRIDobj of elevations.
%   A: GRIDobj of cumulative drainage areas.
%   DB: GRIDobj of drainage basins.
%   FD: FLOWobj of flow directions.
%   topN: Largest number of basins to analyze.
%   channelThreshold: Drainage area threshold for drainage areas.
%
% Output:
%   Ss: Array of streams associated with the largest topN basins.
%   S_Stats: Concavity statistics of the largest topN basin channels,
%       including drainage area, gradient, ks, and concavity (theta).
%   basinIDs: The IDs of the largest topN basins used for the analysis.

Ss = [];
S_Stats = [];

keep = zeros(size(bIAL(:,1)));
for i = 1:size(bIAL,1)
    if sum(basinIDs == bIAL(i,1)) > 0
        keep(i) = 1;
    end
end

bIAL(keep==0,:) = [];
bIAL = sortrows(bIAL,2,'descend');


% bIAL = sortrows(bIAL,2,'descend');
% bIAL(topN+1:end,:) = [];
% basinIDs = bIAL(:,1);

[dbG,~,~] = GRIDobj2mat(DB);

for i = 1:size(bIAL,1)
    dbT = dbG==bIAL(i,1);
    [FDc,~] = crop(FD,dbT);
    [DEMc,~] = crop(DEM,dbT);
    [Ac,~] = crop(A,dbT);
    
    try
        S = STREAMobj(FDc,'minarea',channelThreshold,'unit','mapunits');
        Stats = slopearea(S,DEMc,Ac,'plot',false);
        if abs(Stats.theta) > 2
            Stats.a = NaN;
            Stats.g = NaN;
            Stats.ks = NaN;
            Stats.theta = NaN;
        end
    catch er
        Stats.a = NaN;
        Stats.g = NaN;
        Stats.ks = NaN;
        Stats.theta = NaN;
    end
    
    Ss = [Ss;S];
    S_Stats = [S_Stats;Stats];
end