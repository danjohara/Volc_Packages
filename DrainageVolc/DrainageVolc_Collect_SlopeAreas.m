function [allAS,allFlowAS,DBxy,DBxy_cell,DB_ID] = DrainageVolc_Collect_SlopeAreas(bIAL,DB,D,FD,S,A,dx,topN,DEM)
% Name: CollectSlopeAreas
% Author: Daniel O'Hara
% Date: 05/11/2021 (mm/dd/yyyy)
% Description: Script to collect the contributing drainage area and slopes
%   of the largest basins.
%
% Input:
%   bIAL: Array of basin statistics from Collect_DrainageBasin_Stats_wCross.
%   DB: GRIDobj of drainage basins.
%   D: GRIDobj of basin flow distances.
%   FD: FLOWobj of flow directions.
%   S: GRIDobj of slopes.
%   A: GRIDobj of cumulative drainage areas.
%   dx: Grid resolution.
%   topN: Largest number of basins to analyze (see comments on 
%       DrainageVolc_Analysis for values.
%   DEM: GRIDobj of topography.
%
% Output:
%   allAS: Array that contains all of the drainage area and slope values of
%       the largest topN basins.
%   allFlowAS: Array that contains the drainage area and slope values of
%       only the longest flow paths within the largest topN basins.
%   DBxy: N x 2 array of largest drainage basin X- and Y- coordinates, 
%       separated by NaN's for each basin.
%   DBxy_cell: Cell array of largest drainage basin X- and Y- coordinates.
%   DB_ID: Array of largest basin IDs.

allAS = [];
allFlowAS = [];

[dbG,x,y] = GRIDobj2mat(DB);
[sG,~,~] = GRIDobj2mat(S);
[aG,~,~] = GRIDobj2mat(A);
[dG,~,~] = GRIDobj2mat(D);
[zG,~,~] = GRIDobj2mat(DEM);

if topN >= 1
    bIAL = sortrows(bIAL,2,'descend');
    bIAL(topN+1:end,:) = [];
elseif topN > 0
    bIAL = sortrows(bIAL,2,'descend');

    useNum = ceil(size(bIAL,1)*(topN));
    bIAL(useNum+1:end,:) = [];
elseif topN < 0 && topN > -1
    tmpDBG = dbG;
    zG = zG - min(zG(:));
    zG = zG./max(zG(:));

    tmpDBG(zG < (1-abs(topN))) = NaN;
    
    tmpIDs = unique(tmpDBG(:));
    
    keep = zeros(size(bIAL(:,1)));
    for i = 1:size(bIAL,1)
        if sum(tmpIDs == bIAL(i,1)) > 0
            keep(i) = 1;
        end
    end
    
    bIAL(keep==0,:) = [];
    bIAL = sortrows(bIAL,2,'descend');
else
    error('TopN has incorrect values.')
end

DB_ID = bIAL(:,1);
DBxy = [];
DBxy_cell = {};

for i = 1:size(bIAL,1)
    dbT = dbG==bIAL(i,1);
    
    sT = sG;
    sT(dbT==0) = NaN;
    
    aT = aG.*dx^2;
    aT(dbT==0) = NaN;
    
    asT = [aT(:),sT(:)];
    tt = isnan(sum(asT,2));
    asT(tt,:) = [];
    allAS = [allAS;asT(:,1),asT(:,2)];
    
    tmpD= dG;
    tmpD(dbT==0) = NaN;
    [maxDi,maxDj] = find(tmpD==max(tmpD(:)),1);
    
    maxDLin = sub2ind(size(tmpD),maxDi,maxDj);
    [channelI,~,~,~] = flowpathextract(FD,maxDLin);
    
    for j = 1:length(channelI)
        allFlowAS = [allFlowAS;aT(channelI(j)),sT(channelI(j))];
    end
    
    bb = bwboundaries(dbT);
    tmp = [];
    for j= 1:size(bb{1},1)
        tmp = [tmp;x(bb{1}(j,2)),y(bb{1}(j,1))];
    end
    DBxy = [DBxy;tmp;NaN,NaN];
    DBxy_cell = [DBxy_cell;{tmp}];
end