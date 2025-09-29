function [ei,ii,contDi,elDi] = MorVolc_EI_II(contXY,el,closedContour,bfe)
% Name: MorVolc_Ellipticity_Irregularity
% Author: Daniel O'Hara
% Date: 02/24/2021 (mm/dd/yyyy)
% Description: Script to calculate the ellipticity and irregularity index 
%   of a best-fitting ellipse defined by a set of points.
%
% Input:
%   contXY: Array or cell array of contour xy values. If closedContour is 
%       set to 1, contXY is assumed to be a Nx2 array of x-y values. If
%       closedContour is set to 0, contXY is assumed to be a cell array of 
%       multiple entries.
%   el: Best-fitting ellipse.
%   closedContour: Flag to indicate if this contour is closed or partial.
%   bfe: Flag to indicate whether this is using maximum diameter (like the
%       original MORVOLC) or a best-fitting ellipse.
%
% Output:
%   ei: Ellipticity index.
%   ii: Irregularity index.
%   contDi: Boundary dissection index.
%   elDi: Best-fitting ellipse dissection index.

%% Closed Contour algorithm
if closedContour
    %% Get Points
    x = contXY(:,1);
    y = contXY(:,2);
    
    %% Ellipticity index
    [Area,~] = MorVolc_Calculate_Area_Width(x,y);
    if bfe
        L = el.longAxis;
    else
        pd = pdist2([x,y],[x,y]);
        L = max(pd(:))/2;
    end

    ei = (pi*L^2)/Area;

    %% Contour dissection index
    contPerim = 0;
    xy = [x,y];
    for j = 1:length(x)-1
        contPerim = contPerim + norm(xy(j,:) - xy(j+1,:));
    end
    contPerim = contPerim + norm(xy(end,:)-xy(1,:));
    
    contDi = contPerim/(2*Area)*sqrt(Area/pi);
    
    %% Best-fitting ellipse dissection index.
    A = pi*L^2/ei;
    b = A/pi/L;
    elPerim = pi*(3*(L+b)-sqrt((3*L+b)*(L+3*b)));
    
    elDi = elPerim/(2*Area)*sqrt(Area/pi);
    
    %% Irregularity Index
    ii = contDi-(elDi-1);
else
    acc = .001;
    azConnectThresh = 5;

    %% Get Points and Azimuths
    x = 0;
    y = 0;
    az = NaN;

    if ~iscell(contXY)
        contXY = {contXY};
    end
    
    for i= 1:length(contXY)
        x = [x;contXY{i}(:,1)-el.x0;0];
        y = [y;contXY{i}(:,2)-el.y0;0];

        tmpAz = [];
        for j = 1:size(contXY{i},1)
            phi = Calculate_Azimuths([el.x0,el.y0],contXY{i}(j,1:2));
            tmpAz = [tmpAz;phi];
        end

        az = [az;tmpAz;NaN];
    end

    %% Correct points to help remove odd geometries
    ii = find(x == 0);
    for i = length(ii)-1:-1:2
        if abs(az(ii(i)+1) - az(ii(i)-1)) <= azConnectThresh
            x(ii(i)) = [];
            y(ii(i)) = [];
            az(ii(i)) = [];
        end
    end

    %% Calculate azimuth ranges.
    ii = find(x==0);
    azRanges = [];
    for i = 2:length(ii)
        tmpAzs = az(ii(i-1)+1:ii(i)-1);
        azRanges = [azRanges;min(tmpAzs),max(tmpAzs)];
    end
    azRanges = sortrows(azRanges,1);

    %% Ellipticity index
    Area = polyarea(x,y);
    if bfe
        L = el.longAxis;
    else
        pd = sqrt(x.^2 + y.^2);
        L = max(pd(:));
    end

    circX = 0;
    circY = 0;
    for i = 1:size(azRanges,1)
        useDegs = azRanges(i,1):acc:azRanges(i,2);
        tmpX = L*sin(pi*2*(useDegs)/360);
        tmpY = L*cos(pi*2*(useDegs)/360);

        circX = [circX;tmpX';0];
        circY = [circY;tmpY';0];
    end

    circA = polyarea(circX,circY);
    ei = circA/Area;

    %% Contour dissection index
    contPerim = 0;
    xy = [x,y];
    for j = 1:length(x)-1
        contPerim = contPerim + norm(xy(j,:) - xy(j+1,:));
    end
    contPerim = contPerim + norm(xy(end,:)-xy(1,:));
    
    contDi = contPerim/(2*Area)*sqrt(Area/pi);
    
    %% Best-fitting ellipse dissection index.
    A = pi*L^2/ei;
    b = A/pi/L;
    elPerim = pi*(3*(L+b)-sqrt((3*L+b)*(L+3*b)));
    
    elDi = elPerim/(2*Area)*sqrt(Area/pi);
    
    % circPerim = 0;
    % cxy = [circX,circY];
    % for j = 1:length(x)-1
    %     circPerim = circPerim + norm(cxy(j,:) - cxy(j+1,:));
    % end
    % circPerim = circPerim + norm(cxy(end,:)-cxy(1,:));
    % 
    % elDi = circPerim/(2*circA)*sqrt(circA/pi);

    %% Irregularity Index
    ii = contDi-(elDi-1);
end
end
