function [surfXYs,surfZs,maxDepth,maxVol,volcDepths,volcVols] = Morvolc_CalculateHeightVols_Crater(X,Y,Z,craterXYZ,exZ,interpSurfaces)
% Name: Morvolc_CalculateHeightVols_Crater
% Author: Daniel O'Hara
% Date: 03/17/2021 (mm/dd/yyyy)
% Description: Script to calculate the volume, depth, maximum volume, and
%   maximum depth of a given crater.
%
% Input:
%   X: Grid of X-coordinates.
%   Y: Grid of Y-coordinates.
%   Z: Grid of elevations.
%   craterXYZ: XYZ coordinates of the crater boundary.
%   exZ: Extreme value for maximum volume/height (maximum crater boundary 
%       elevation).
%   interpSurfaces: Class of logical values to determine interpolation
%       method of the crater surface, includes:
%           Natural: Matlab natural interpolation.
%           IDW: Inverse Distance Weighting interpolation.
%           Kringing: Kringing interpolation.
%           GreensFunc: Interpolation using a Green's function.
%
% Output:
%   surfXYs: X-Y interpolation points for basal surface.
%   surfZs: Interpolated Z values for the surfXYs value pairs, given as a
%       class which follows the same structure as interpSurfaces.
%   maxDepth: Maximum depth, calculated as the distance between the 
%       crater minima and the exZ elevation.
%   maxVol: Maximum volume, calculated as the volume between the crater 
%       maxima and the exZ elevation.
%   volcDepths: Crater depths, calculated as the distance between the 
%       crater maxima and interpolated crater surface of the crater 
%       boundary, given as a class which follows the same structure as 
%       interpSurfaces.
%   volcVols: Edifice volume, calculated as the amount of material between 
%       the crater maxima and interpolated crater surface of the crater 
%       boundary, given as a class which follows the same structure as 
%       interpSurfaces.

%% Find Maximum Height & Location
[ii,jj] = find(Z==min(Z(:)),1);
lowHeight = Z(ii,jj);
lowX = X(ii,jj);
lowY = Y(ii,jj);

dx = X(2,2)-X(1,1);

%% Determine Grid Points in Boundary
xx = X(:);
yy = Y(:);
zz = Z(:);

[inp,~] = inpolygon(xx,yy,craterXYZ(:,1),craterXYZ(:,2));
surfXYs = [xx(inp),yy(inp)];

xxp = xx(inp);
yyp = yy(inp);
t1 = xxp==lowX;
t2 = yyp==lowY;
t3 = t1.*t2;
peakI = find(t3==1,1);

%% Generate Surfaces & Get Height/Volumes
if interpSurfaces.Natural
    boundaryInterp = scatteredInterpolant(craterXYZ(:,1),craterXYZ(:,2),craterXYZ(:,3),'natural');
    basalSurfZ = boundaryInterp(xxp,yyp);
    
    surfZs.Natural = basalSurfZ;
    volcDepths.Natural = abs(lowHeight-basalSurfZ(peakI));
    
    volcVols.Natural = sum(abs(basalSurfZ-zz(inp)))*dx^2;
else
    surfZs.Natural = [];
    volcDepths.Natural = [];
    volcVols.Natural = [];
end

if interpSurfaces.IDW
    basalSurfZ = idw(craterXYZ(:,1:2),craterXYZ(:,3),[xxp,yyp]);
    
    surfZs.IDW = basalSurfZ;
    volcDepths.IDW = abs(lowHeight-basalSurfZ(peakI));
    
    volcVols.IDW = sum(abs(basalSurfZ-zz(inp)))*dx^2;
else
    surfZs.IDW = [];
    volcDepths.IDW = [];
    volcVols.IDW = [];
end

if interpSurfaces.Kringing
    v = variogram(craterXYZ(:,1:2),craterXYZ(:,3),'plotit',false);
    [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable','plotit',false);
    [basalSurfZ,~] = kriging(vstruct,craterXYZ(:,1),craterXYZ(:,2),craterXYZ(:,3),xxp,yyp);

    surfZs.Kringing = basalSurfZ;
    volcDepths.Kringing = abs(lowHeight-basalSurfZ(peakI));
    
    volcVols.Kringing = sum(abs(basalSurfZ-zz(inp)))*dx^2;
else
    surfZs.Kringing = [];
    volcDepths.Kringing = [];
    volcVols.Kringing = [];
end

%% Calculate Max Height/Vol
maxDepth = abs(min(Z(:))-exZ);

basalSurfZ = ones(size(xx(inp)))*exZ;

maxVol = sum(abs(basalSurfZ-zz(inp)))*dx^2;