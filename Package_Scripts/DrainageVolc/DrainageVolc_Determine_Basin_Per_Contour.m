function Basin_Contour_ContourP_Count_Length_Area = DrainageVolc_Determine_Basin_Per_Contour(DEMf,DB,A,basinContIter,areaThreshold)
% Name: Determine_Basin_Per_Contour
% Date: 05/19/2022 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to determine the number of baisns that exists at
% specific contours of the edifice.
%
% Input:
%   DEMf: DEM clipped to the edifice boundary (GRIDobj).
%   DB: GRIDobj of drainage basins.
%   A: GRIDobj of maximum basin drainage areas.
%   basinContIter: Contour interval for analysis, can be either in meters,
%       or as a percent (given as a negative value between -1 and 0).
%   areaThreshold: Drainage area threshold for secondary basin count.
%
% Output:
%   Basin_Contour_ContourP_Count_Length_Area: Array containing analyzed
%       contour values (in meters and percent relief from main flank), the 
%       number of basins in each contour, the length of the contour, and 
%       the area of each contour.

%% Setup
Basin_Contour_ContourP_Count_Length_Area = [];

[Z,x,y] = GRIDobj2mat(DEMf);
if size(x,1) == 1
    x = x';
end

[DBg,~,~] = GRIDobj2mat(DB);

DBg(isnan(Z)) = NaN;

if ~isnan(areaThreshold)
    [Ag,~,~] = GRIDobj2mat(A);
    DBg(Ag<areaThreshold) = NaN;
end

tmpZ = double(Z);

%% Scale edifice
tmpMinZ = min(tmpZ(:));
tmpZ = tmpZ - tmpMinZ;
tmpMaxZ = max(tmpZ(:));

if basinContIter < 0 && basinContIter > -1
    isPer = 1;
    tmpZ = tmpZ./tmpMaxZ;
else
    isPer = 0;
end

% tmpZ(isnan(tmpZ)) = -1;
%% Determine number of basins per contour
basinContIter = abs(basinContIter);

for i = 0:basinContIter:max(tmpZ(:))
    if isPer
        contVal_Per = i;
        contVal = tmpMaxZ*i + tmpMinZ;
    else
        contVal = i + tmpMinZ;
        contVal_Per = i/tmpMaxZ;
    end

    tmpDB = DBg;
    tmpDB(tmpZ<i) = NaN;
    uniBasins = unique(tmpDB(:));
    uniBasins(isnan(uniBasins)) = [];
    basinCount = numel(uniBasins);
    contArea = sum(~isnan(tmpDB(:)))*DEMf.cellsize^2;
    
    try
        tmpZ2 = tmpZ;
        tmpZ2(tmpZ<i) = NaN;

        bb = bwboundaries(~isnan(tmpZ2));
        contLength = 0;
        for j = 1:length(bb)
            cx = x(bb{j}(:,2));
            cy = y(bb{j}(:,1));
            cxy = [cx,cy];

            diffCC = [diff(cxy,1);cxy(end,:)-cxy(1,:)];
            contLength = contLength + sum(sqrt(sum(diffCC .* diffCC, 2)));
        end
    catch er
        contLength = NaN;
    end

    Basin_Contour_ContourP_Count_Length_Area = [Basin_Contour_ContourP_Count_Length_Area;
        contVal,contVal_Per,basinCount,contLength,contArea];
end
