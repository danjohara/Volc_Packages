function [surfXYs,surfZs,maxHeight,maxVol,volcHeights,volcVols,anyHeights] = MorVolc_CalculateHeightVols(X,Y,Z,boundaryXYZ,exZ,interpSurfaces)
% Name: MorVolc_CalculateHeightVols
% Author: Daniel O'Hara
% Date: 03/17/2021 (mm/dd/yyyy)
% Description: Script to calculate the volume, height, maximum volume, and
%   maximum height of a given edifice.
%
% Input:
%   X: Grid of X-coordinates.
%   Y: Grid of Y-coordinates.
%   Z: Grid of elevations.
%   boundaryXYZ: XYZ coordinates of the edifice boundary.
%   exZ: Extreme value for maximum volume/height (minimum edifice boundary 
%       elevation).
%   interpSurfaces: Class of logical values to determine interpolation
%       method of the basal surface, includes:
%           Natural: Matlab natural interpolation.
%           IDW: Inverse Distance Weighting interpolation.
%           Kriging: Kriging interpolation.
%
% Output:
%   surfXYs: X-Y interpolation points for basal surface.
%   surfZs: Interpolated Z values for the surfXYs value pairs, given as a
%       class which follows the same structure as interpSurfaces.
%   maxHeight: Maximum height, calculated as the distance between the 
%       edifice peak and the exZ elevation.
%   maxVol: Maximum volume, calculated as the volume between the edifice 
%       peak and the exZ elevation.
%   volcHeights: Edifice height, calculated as the distance between the 
%       edifice peak and interpolated basal surface of the edifice 
%       boundary, given as a class which follows the same structure as 
%       interpSurfaces.
%   volcVols: Edifice volume, calculated as the amount of material between the 
%       edifice peak and interpolated basal surface of the edifice 
%       boundary, given as a class which follows the same structure as 
%       interpSurfaces.

%% Find Maximum Height & Location
[ii,jj] = find(Z==max(Z(:)),1);
peakHeight = Z(ii,jj);
peakX = X(ii,jj);
peakY = Y(ii,jj);

dx = X(2,2)-X(1,1);
anyHeights = [];

%% Fix Boundary if needed
ff = find(isnan(boundaryXYZ(:,1)),1);
if ~isempty(ff)
    boundaryXYZ(ff:end,:) = [];
end

boundaryXYZ(isnan(boundaryXYZ(:,3)),:) = [];

%% Determine Grid Points in Boundary
xx = X(:);
yy = Y(:);
zz = Z(:);

[inp,~] = inpolygon(xx,yy,boundaryXYZ(:,1),boundaryXYZ(:,2));
surfXYs = [xx(inp),yy(inp)];

xxp = xx(inp);
yyp = yy(inp);
zzp = zz(inp);
t1 = xxp==peakX;
t2 = yyp==peakY;
t3 = t1.*t2;
peakI = find(t3==1,1);

%% Generate Surfaces & Get Height/Volumes
if interpSurfaces.Natural
    boundaryInterp = scatteredInterpolant(boundaryXYZ(:,1),boundaryXYZ(:,2),boundaryXYZ(:,3),'natural');
    basalSurfZ = boundaryInterp(xxp,yyp);
    
    surfZs.Natural = basalSurfZ;
    volcHeights.Natural = peakHeight-basalSurfZ(peakI);
    
    tmp = zzp-basalSurfZ;
    [ii] = find(tmp==max(tmp(:)),1);
    anyHeights.Natural.X = xxp(ii);
    anyHeights.Natural.Y = yyp(ii);
    anyHeights.Natural.Height = tmp(ii);
    
    volcVols.Natural = sum(abs(basalSurfZ-zzp))*dx^2;
else
    surfZs.Natural = [];
    volcHeights.Natural = [];
    volcVols.Natural = [];
    anyHeights.Natural = [];
end

if interpSurfaces.IDW
    basalSurfZ = idw(boundaryXYZ(:,1:2),boundaryXYZ(:,3),[xxp,yyp]);
    
    surfZs.IDW = basalSurfZ;
    volcHeights.IDW = peakHeight-basalSurfZ(peakI);
    
    tmp = zzp-basalSurfZ;
    [ii] = find(tmp==max(tmp(:)),1);
    anyHeights.IDW.X = xxp(ii);
    anyHeights.IDW.Y = yyp(ii);
    anyHeights.IDW.Height = tmp(ii);
    
    volcVols.IDW = sum(abs(basalSurfZ-zzp))*dx^2;
else
    surfZs.IDW = [];
    volcHeights.IDW = [];
    volcVols.IDW = [];
    anyHeights.IDW = [];
end

if interpSurfaces.Kriging
    v = variogram(boundaryXYZ(:,1:2),boundaryXYZ(:,3),'plotit',false);
    [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable','plotit',false);
    [basalSurfZ,~] = kriging(vstruct,boundaryXYZ(:,1),boundaryXYZ(:,2),boundaryXYZ(:,3),xxp,yyp);

    surfZs.Kriging = basalSurfZ;
    volcHeights.Kriging = peakHeight-basalSurfZ(peakI);
    
    tmp = zzp-basalSurfZ;
    [ii] = find(tmp==max(tmp(:)),1);
    anyHeights.Kriging.X = xxp(ii);
    anyHeights.Kriging.Y = yyp(ii);
    anyHeights.Kriging.Height = tmp(ii);
    
    volcVols.Kriging = sum(abs(basalSurfZ-zzp))*dx^2;
else
    surfZs.Kriging = [];
    volcHeights.Kriging = [];
    volcVols.Kriging = [];
    anyHeights.Kriging = [];
end

%% Calculate Max Height/Vol
maxHeight = max(Z(:))-exZ;

basalSurfZ = ones(size(xx(inp)))*exZ;

maxVol = sum(abs(basalSurfZ-zz(inp)))*dx^2;