function [summitXYZ,summitCont] = Morvolc_GetSummitRegion(summitRegion,Z,S,X,Y,lowFlankZ,conts,useMaxZ,peakDiff)
% Name: Morvolc_GetSummitRegion
% Author: Daniel O'Hara
% Date: 03/17/2021 (mm/dd/yyyy)
% Description: Script to determine the summit region of a volcano, defined
%   as either the contour that encompasses multiple peaks of the edifice,
%   or the contour that has the greatest slope change.
%
% Input:
%   summitRegion: Flag for summit region designation.
%   Z: Grid of elevations.
%   S: Grid of slope values.
%   X: Grid of x-coordinates.
%   Y: Grid of y-coordinates.
%   lowFlankZ: Elevation of the lower flank.
%   conts: Array of contours to use for analysis.
%   useMaxZ: Maximum elevation to generate contours for summit analysis.
%   peakDiff: Percentage value to distinguish edifice peaks from
%       local-minima topography.
%
% Output:
%   summitXYZ: Summit boundary x-, y-, and z-coordinates.
%   summitCont: Contour elevation of summit.

%% Setup
x = X(1,:);
y = Y(:,1);
dx = x(2)-x(1);

%% Collect Summit
if summitRegion < 0 && summitRegion > -1
    tmp = nanmin(Z(:)) + range(Z(:))*(1+summitRegion);
    t1 = find(abs(conts-tmp) == min(abs(conts-tmp)),1);
    summitCont = conts(t1);
elseif summitRegion==4 || summitRegion == 5
    %% Check if crater is being used as the summit
    tmp = conts<useMaxZ;
    t1 = find(tmp==1,1,'last');
    summitCont = conts(t1);
elseif summitRegion==3
    %% Determine number of basins per contour length and find highest decrease
    DEM = GRIDobj(X,Y,Z);
    FD = FLOWobj(DEM);
    DB = drainagebasins(FD);

    Basin_Contour_ContourP_Count_Length_Area = DrainageVolc_Determine_Basin_Per_Contour(DEM,DB,-.05);
    Basin_Contour_ContourP_Count_Length_Area(Basin_Contour_ContourP_Count_Length_Area(:,2)>=.9,:) = [];

    contLength = Basin_Contour_ContourP_Count_Length_Area(:,3)./Basin_Contour_ContourP_Count_Length_Area(:,4);
    diffcontLength = [0;diff(contLength)];

    ii = find(diffcontLength<0,1,'last');
    summitCont = Basin_Contour_ContourP_Count_Length_Area(ii,1);

else
    %% Generate mean slope array.
    conts(conts>useMaxZ) = [];
    conts(conts<lowFlankZ) = [];
    conts = [conts,ones(size(conts))];
    meanS_multiBW = zeros(length(conts),2);

    %% Loop through contours 
    for i = 2:length(conts)
        disp(sprintf('   Contour %d / %d',i,length(conts)))
        
        %% Isolate elevations   
        tmpZ = Z;
        tmpZ(tmpZ<conts(i)) = NaN;
        bb = bwboundaries(~isnan(tmpZ),'noholes');
        
        % If multiple groups of topography exist at the contour, check if they 
        %   are actually peaks by comparing their volumes to the main
        %   structure.
        if length(bb) > 1
            rel = zeros(size(bb));
            for j = 1:length(bb)
                px = [];
                py = [];
                for k = 1:size(bb{j},1)
                    px = [px;x(bb{j}(k,2))];
                    py = [py;y(bb{j}(k,1))];
                end
                pp = inpolygon(X,Y,px,py);
                zz = Z;
                zz(~pp) = NaN;
                zz = zz-min(zz(:));
                rel(j) = range(zz(:));
            end
            
            conts(i,2) = find(rel==max(rel),1);
            vC = rel./rel(conts(i,2));
            vC(vC==1) = NaN;
            if sum((vC > peakDiff)*1) > 0
                meanS_multiBW(i,2) = 1;
                break;
            end
        end
        
        % If multiple peaks don't exist, collect the mean slopes of the
        %   contour.
        
        ccI = contourc(Z,[conts(i),conts(i)]);
        ccI = Convert_Contours(ccI,0);
        
        ss = S;
        ss(Z<=conts(i-1)) = NaN;
        ss(Z>conts(i)) = NaN;
        
        meanS_multiBW(i,1) = nanmean(ss(:));
    end
    meanS_multiBW(1,1) = meanS_multiBW(2,1);

    %% Determine summit contour
    % Either find where the contour where multiple peaks exist, or the contour 
    %   that has the highest slope.
    if summitRegion == 2 && sum(meanS_multiBW(:,2)) > 0
        ii = find(meanS_multiBW(:,2)==1);
    else
        tt = gradient(atand(meanS_multiBW(:,1)),conts(2,1)-conts(1,1));
        ii = find(tt==min(tt(2:end)),1);
    end
    
    summitCont = conts(ii,1);
end

%% Collect summit contour points
tmpZ = Z;
tmpZ(tmpZ<summitCont) = NaN;
bb = bwboundaries(~isnan(tmpZ));
summitXYZ = [];

useI = 1;
for i = 2:length(bb)
    if size(bb{i},1) > size(bb{useI},1)
        useI = i;
    end
end

for i = 1:size(bb{useI},1)
    summitXYZ = [summitXYZ;...
        X(bb{useI}(i,1),bb{useI}(i,2)),...
        Y(bb{useI}(i,1),bb{useI}(i,2))];
end
zz = interp2(X,Y,Z,summitXYZ(:,1),summitXYZ(:,2));
summitXYZ = [summitXYZ,zz];
    