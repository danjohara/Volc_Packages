function [cont_area_convHullArea,volumeDiff,convZ] = CalculateErodedVolume(X,Y,Z,bxyz,contZ,startZ,onlyLargestContour,silentRun)
%%
% Name: CalculateErodedVolume
% Author: Daniel O'Hara
% Original Date: 11/2020 (mm/yyyy)
% Updated Date: 03/03/2021 (mm/dd/yyyy)
% Updated Date: 04/18/2023 (mm/dd/yyyy)
% Description: Script to calculate the minimum amount of eroded volume from 
%   an edifice by calculating the areal difference between actual
%   topography and a convex hull, integrated across contours.
%
% Citation: O'Hara, D., Karlstrom, L. (Dec. 2020). Volcanic Edifice Erosion
%   as a Probe of Landscape Evolution in the Cascades Arc. American 
%   Geophysical Union Annual Meeting, Virtual Conference.
%
% Input:
%   X: Grid of x-coordinates.
%   Y: Grid of y-coordinates.
%   Z: Grid of elevations.
%   contZ: Contour interval for analysis.
%   startZ: Starting contour for analysis.
%   onlyLargestContour: Consider the interpolation only by the largest closed
%       contour.
%   silentRun: Flag to run script silently.
%
% Output: 
%   cont_area_convHullArea: Matrix contour elevations, area, and convex
%       hull area.
%   volumeDiff: Eroded volume (volume difference between actual topography
%       and convex hull topography).
%   coneZ: Interpolated surface.

%% Setup
cont_area_convHullArea = [];
counter = 1;
contourInterpolationStep = 20;
gridDx = sqrt((X(1,1)-X(2,2))^2 + (Y(1,1)-Y(2,2))^2);

%% Get contour areas
convPoints = [];
for i = startZ:contZ:max(Z(:))+contZ
    if ~silentRun
        disp(sprintf('%d / %d',counter,length(startZ:contZ:max(Z(:))+contZ)))
    end
    tmpZ = Z>=i;
    
    cutI1 = find(sum(tmpZ,2)>0,1,'first');
    cutI2 = find(sum(tmpZ,2)>0,1,'last');
    cutJ1 = find(sum(tmpZ,1)>0,1,'first');
    cutJ2 = find(sum(tmpZ,1)>0,1,'last');
    
    cutZ = tmpZ(cutI1:cutI2,cutJ1:cutJ2);
    cutX = X(cutI1:cutI2,cutJ1:cutJ2);
    cutY = Y(cutI1:cutI2,cutJ1:cutJ2);
    
    bb = bwboundaries(cutZ);
    if onlyLargestContour
        bbSizes = zeros(size(bb));
        for ii = 1:length(bb)
            bbSizes(ii) = size(bb{ii},1);
        end
        if ~isempty(bb)
            bb = bb(find(bbSizes == max(bbSizes),1));
        end
    end
    
    allArea = 0;
    allConvArea = 0;
    
    for j = 1:length(bb)
        try
            hull = convhull(bb{j}(:,2),bb{j}(:,1));
        catch er
            continue;
        end
        
        xx = zeros(length(bb{j}),1);
        yy = xx;
        for k = 1:length(bb{j})
            xx(k) = cutX(bb{j}(k,1),bb{j}(k,2));
            yy(k) = cutY(bb{j}(k,1),bb{j}(k,2));
        end
        allArea = allArea + polyarea(xx,yy);
        
        %Interpolate between points
        xx = [];
        yy = [];
        for k = 1:length(hull)
            tx = cutX(bb{j}(hull(k),1),bb{j}(hull(k),2));
            ty = cutY(bb{j}(hull(k),1),bb{j}(hull(k),2));

            xx = [xx;tx];
            yy = [yy;ty];

            if k < length(hull)
                otherTx = cutX(bb{j}(hull(k+1),1),bb{j}(hull(k+1),2));
                otherTy = cutY(bb{j}(hull(k+1),1),bb{j}(hull(k+1),2));
            else
                otherTx = cutX(bb{j}(hull(1),1),bb{j}(hull(1),2));
                otherTy = cutY(bb{j}(hull(1),1),bb{j}(hull(1),2));
            end
            
            patho = sqrt((tx-otherTx).^2 + (ty-otherTy).^2);
            if patho > contourInterpolationStep
                ddx = tx-otherTx;
                ddy = ty-otherTy;
                slp = ddy/ddx;

                for ii = contourInterpolationStep:contourInterpolationStep:patho
                    newX = tx - sign(ddx)*ii*cos(atan(abs(slp)));
                    newY = ty - sign(ddy)*ii*sin(atan(abs(slp)));
                    xx = [xx;newX];
                    yy = [yy;newY];
                end
%                 disp('here')
            end
        end

        allConvArea = allConvArea + polyarea(xx,yy);
        zz = ones(size(xx))*i;
        convPoints = [convPoints;xx,yy,zz];
    end
    
    cont_area_convHullArea = [cont_area_convHullArea;i,allArea,allConvArea];
    counter = counter + 1;
end

try
    % Remove points outside of the boundary
    pp = inpolygon(convPoints(:,1),convPoints(:,2),bxyz(:,1),bxyz(:,2));
    convPoints(pp==0,:) = [];

%     % Remove duplicate XY points
%     convPoints = sortrows(round(convPoints,0),3,'descend');
%     cutPoints = zeros(size(convPoints,1),1);
%     for i = 1:size(convPoints,1)
%         tmp = pdist2(convPoints(i,1:2),convPoints(:,1:2));
%         t3 = tmp < abs(X(2,2)-X(1,1));
%         ii = find(t3==1);
%         ii(ii<i) = [];
%         if length(ii) > 1
%             cutPoints(ii(2:end)) = 1;
%         end
%     end
%     convPoints(cutPoints==1,:) = [];

    cInterp = scatteredInterpolant(convPoints(:,1),convPoints(:,2),convPoints(:,3),'linear','none');
catch er 
    disp('  TOO FEW POINTS IN GRID')
    volumeDiff = NaN;
    cont_area_convHullArea = NaN;
    convZ = NaN;
    return;
end
convZ = cInterp(X,Y);

bb = bwboundaries(~isnan(Z));
bx = [];
by = [];
for i = 1:size(bb{1},1)
    bx = [bx;X(bb{1}(i,1),bb{1}(i,2))];
    by = [by;Y(bb{1}(i,1),bb{1}(i,2))];
end

inP = inpolygon(X,Y,bx,by);
convZ(~inP) = NaN;

%% Calculate Volume
gridDiff = convZ-Z;
gridDiff(gridDiff<0) = 0;
gridDiff(Z<startZ) = 0;
volumeDiff = nansum(gridDiff(:))*gridDx^2;

end