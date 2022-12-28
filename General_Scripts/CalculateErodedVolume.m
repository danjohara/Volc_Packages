function [cont_area_convHullArea,volumeDiff] = CalculateErodedVolume(X,Y,Z,contZ,startZ,silentRun)
%%
% Name: CalculateErodedVolume
% Author: Daniel O'Hara
% Original Date: 11/2020 (mm/yyyy)
% Updated Data: 03/03/2021 (mm/dd/yyyy)
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
%   silentRun: Flag to run script silently.
%
% Output: 
%   cont_area_convHullArea: Matrix contour elevations, area, and convex
%       hull area.
%   volumeDiff: Eroded volume (volume difference between actual topography
%       and convex hull topography).
%% Setup
cont_area_convHullArea = [];
counter = 1;

%% Get contour areas
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
        
        xx = zeros(length(hull),1);
        yy = xx;
        for k = 1:length(hull)
            xx(k) = cutX(bb{j}(hull(k),1),bb{j}(hull(k),2));
            yy(k) = cutY(bb{j}(hull(k),1),bb{j}(hull(k),2));
        end
        allConvArea = allConvArea + polyarea(xx,yy);
    end
    
    cont_area_convHullArea = [cont_area_convHullArea;i,allArea,allConvArea];
    counter = counter + 1;
end

%% Calculate Volume
trueVol = trapz(cont_area_convHullArea(:,2))*contZ;
convHullVol = trapz(cont_area_convHullArea(:,3))*contZ;
volumeDiff = convHullVol-trueVol;
end