function morRes = MorVolc_Analysis(pack)
%%
% Name: MorVolc_Analysis
% Author: Daniel O'Hara
% Original MorVolc Author: Pablo Grosse
% Date: 02/24/2021 (mm/dd/yyyy)
% Description: Script to calculate the MorVolc statistics of an edifice 
%   (Grosse et al., 2009, 2012).  
% 
% NOTE:This script uses TopoToolbox (https://github.com/wschwanghart/topotoolbox), 
%   as well as some other community-created functions (e.g., fitellipse, 
%   ll2utm) and requires Matlab's Mapping, Parallel Computing, and Image 
%   Processing Toolboxes.
%
% ALSO: This script is continuously being updated. Some functionality is 
%   implemented in MorVolc that needs fixed (e.g., peak counts); this is 
%   underway. Users are free to use and adapt this script as needed.
%
% Input:
%   pack: Strucure of input parameters, includes:
%       tifFile: Tif file and path of edifice DEM.
%       boundaryXY: Array, or .shp file and path, of edifice boundary X- and 
%           Y-coordinates.
%       maskMap: Cell array, or .shp file and path, of mask regions X- and 
%           Y-coordinates to ignore in analysis.
%       craterXY: Array, or .shp file and path, of edifice crater X- and 
%           Y-coordinates.
%
%       dx: Resolution of DEM (in meters).
%       peakDiff: Percentage volume difference to distinguish edifice peaks
%           from local topographic maxima.
%       contIter: Numerical contour interval value for contour analysis.
%           If given as a negative value between -1 and 0, value will be
%           used as a percentage of the total edifice height. If contIter 
%           is an array, the algorithm will use the given values as the
%           contours (negative flag still applies).
%       craterContIter: Numerical contour interval value for crater contour
%           analysis. If given as a negative value between -1 and 0, value 
%           will be used as a percentage of the total crater height.
%       peakContIter: Numerical contour interval value for peak count
%           analysis. If given as a negative value between -1 and 0, value 
%           will be used as a percentage of the total edifice height
%
%       correctCrater: Flag to attempt to correct the crater region and pin
%           it to the surrounding divides (requires a closed crater).
%       summitRegion: Flag to determine the algorithm used to define the
%           summit region:
%               1: Summit defined by the greatest slope change.
%               2: Summit defined by the contour that encompasses multiple
%                   peaks.
%               3: Summit defined by the upper-most decrease in number of
%                   drainage basins per contour (ignoring upper 10% of the
%                   edifice).
%               4: Summit defined by the crater.
%               Value raging -1 to 0: Summit defined by upper percentile of
%                   edifice height (similar to DrainageVolc basinTopN flag.
%       ignoreSummitIndex: Flag to ignore summit and perform ellipticity
%           and irregularity analysis over edifice.
%
%       interpSurfaces: Structure to calculate the basal and crater
%           surfaces using a variety of interpolation values, including:
%               Natural: Matlab natural interpolation.
%               IDW: Inverse Distance Weighting interpolation.
%               Kringing: Kringing interpolation.
%
%       xlsFile: MS Excel file to save analysis results. Leave empty to not
%           save values.
%       saveFigFolder: Folder to save analysis plots (must include '/' or 
%           '\' at the end, and plotResults must be set to 1 to use). Set 
%           to empty if not saving plots.
%
%       plotResults: Flag to plot analysis results.
%       figPrefix: Name prefix for all figures.
%       saveResFolder: Folder to save output structure. Set to empty if not
%           saving results.
%
% Output:
%   morRes: Analysis results given as a structure. Fields are:
%       GeneralParams:
%           inputs: Script input structure.
%           X: Grid of DEM x-coordinates.
%           Y: Grid of DEM y-coordinates.
%           Z: Grid of DEM elevations.
%           DEM: GRIDobj of DEM (used by TopoToolbox).
%           S: Grid of DEM slope values.
%           C: Grid of DEM profile curvature values.
%           BoundaryXYZ: Matrix of edifice boundary x, y, and z values.
%           Summit_Contour: Elevation between main flank and summit.
%           SummitXYZ: Matrix of summit boundary x, y, and z values.
%           CraterXYZ: Matrix of crater boundary x, y, and z values.
%           MaskXY: Mask region x- and y-coordinates.
%           Basal_Surface_XY: X-Y coordinates of the edifice's basal
%               surface.
%           Basal_Surface_Z: Basal surface interpolated elevations for the
%               X-Y coordinates in Basal_Surface_XY, given as a structure 
%               which follows the interpSurfaces input.
%           Lower_Flank_Contour: Elevation between main flank and lower
%               flank.
%           Lower_Flank_XY: X-Y coordinates of the lower flank boundary.
%           Version: Analysis script version number.
%           StartTime: DateTime object giving the time the script was
%               started.
%           EndTime: DateTime object giving the time the script ended.
%       SizeParams:
%           Basal_Area: Area of edifice boundary.
%           Basal_Width: Diameter of circle with area equivilant to
%               edifice boundary.
%           Major_Basal_Axis: Best-fitting ellipse major axis length of
%               edifice boundary.
%           Minor_Basal_Axis: Best-fitting ellipse minor axis length of
%               edifice boundary.
%           Basal_Axis_Ellipticity: Ellipticity (minor / major axis) value
%               of the edifice boundary.
%           Basel_Ellipse: Best-fitting ellipse of edifice boundary.
%           Summit_Area: Area of summit boundary.
%           Summit_Width: Diameter of circle with area equivilant to
%               summit boundary.
%           Summit_Basal_Axis: Best-fitting ellipse major axis length of
%               summit boundary.
%           Basal_Axis_Ellipticity: Ellipticity (minor / major axis) value
%               of the summit region.
%           Summit_Ellipse:Best-fitting ellipse of summit boundary.
%           Peak_Height: Heights of edifice, defined as distance between 
%               peak and basal surface under peak, given as a structure 
%               which follows the interpSurfaces input.
%           Any_Height: Edifice height anywhere, defined as maximum 
%               distance between edifice and basal surface, given as a 
%               structure which follows the interpSurfaces input.
%           Maximum Height: Maximum height of edifice, defined as distance
%               between peak and lowest boundary elevation.
%           Volume: Volumes of edifice between peak and basel surface, 
%               given as a structure which follows the interpSurfaces 
%               input.
%           Maximum_Volume: Maximum volume of edifice between peak and
%               lowest boundary elevation.
%           Minimum_Eroded_Volume: Minimum edifice eroded volume calculated
%               as the integrated areal difference between actual
%               topography and topography fit by a convex hull over
%               contours.
%           Convex_Hull_Areas: Contour and convex hull contour areas from
%               the eroded volume analysis.
%       OrientParams:
%           Basal_Major_Axis_Azimuth: Azimuthal direction of best-fitting
%               edifice boundary ellipse.
%           Summit_Major_Axis_Azimuth: Azimuthal direction of best-fitting
%               edifice summit ellipse.
%           Contour_Values: Contour values used for contour analysis.
%           Contour_Major_Axis_Azimuth: Azimuthal direction of best-fitting
%               ellipses along each contour.
%       ShapeParams:
%           Contour_Values: Contour values used for contour analysis.
%           Contour_MaxDiam_Ellipticity_Indexes: Ellipticity indexes along 
%               each contour using the maximum diameter of the contour.
%           Contour_BFEllipse_Ellipticity_Indexes: Ellipticity indexes  
%               along each contour using the contour best-fitting ellipse 
%               major axis.
%           Contour_MaxDiam_Irregularity_Indexes: Irregularity indexes 
%               along each contour using the maximum diameter of the 
%               contour.
%           Contour_BFEllipse_Irregularity_Indexes: Irregularity indexes 
%               along each contour using the contour best-fitting ellipse 
%               major axis.
%           Contour_Axis_Ellipticity: Axis ellipticity indexes
%               (best-fitting ellipse minor axis / major axis).
%           Mean_MaxDiam_Ellipticity_Index: Mean contour ellipticity index 
%               using the maximum diameter of the contour.
%           Mean_BFEllipse_Ellipticity_Index: Mean contour ellipticity  
%               index using the contour best-fitting ellipse major axis.
%           Mean_MaxDiam_Irregularity_Index: Mean contour irreguarity index  
%               using the maximum diameter of the contour.
%           Mean_BFEllipse_Irregularity_Index: Mean contour irreguarity 
%               index using the contour best-fitting ellipse major axis.
%           Mean_Ellipse_Ellipticity: Mean ellipse axis ellipticity value.
%           Height_BasalWidth: Ratio between edifice height and basal
%               width.
%           SummitWidth_BasalWidth: Ratio between summit width and basal
%               width.
%           Contour_Ellipses: Best-fitting ellipses along each contour.
%           Skewness: Edifice X-Y skewness.
%           Kurtosis: Edifice X-Y kurtosis.
%       SlopeParams:
%           WholeEdifice_Mean_Slope: Mean slope of entire edifice.
%           WholeEdifice_Median_Slope: Median slope of entire edifice.
%           WholeEdifice_Std_Slope: Standard deviation slope of entire
%               edifice.
%           Flank_Mean_Slope: Mean slope of main edifice flank region.
%           Flank_Median_Slope: Median slope of main edifice flank 
%               region.
%           Flank_Std_Slope: Standard deviation slope of main edifice flank
%               region.
%           Lower_Flank_Mean_Slope: Mean slope of lower edifice flank 
%               region.
%           Lower_Flank_Median_Slope: Median slope of lower edifice flank 
%               region.
%           Lower_Flank_Std_Slope: Standard deviation slope of lower
%               edifice flank region.
%           Summit_Mean_Slope: Mean slope of the summit region.
%           Summit_Median_Slope: Median slope of the summit region.
%           Summit_Std_Slope: Standard deviation slope of summit region.
%           Contour_Values: Contour values used for contour analysis.
%           Contour_Mean_Slopes: Mean slope along each contour.
%           Contour_Median_Slopes: Median slope along each contour.
%           Contour_Min_Slopes: Min slope along each contour.
%           Contour_Max_Slopes: Max slope along each contour.
%           Contour_Max_Mean_Slope: Maximum mean slope of all contours.
%           Contour_Max_Median_Slope: Maximum median slope of all contours.
%           Max_Slope_Height: Elevation of the highest mean slope contour.
%           Max_Slope_Height_Fraction: Percentage of maximum mean slope 
%               elevation compared to entire edifice relief. 
%           Max_Slope: Maximum mean contour slope value.
%       CraterParams:
%           Crater_Area: Area of crater boundary.
%           Crater_Width: Diameter of circle with area equivilant to
%               crater boundary.
%           Crater_Major_Axis: Best-fitting ellipse major axis length of
%               crater boundary.
%           Crater_Minor_Axis: Best-fitting ellipse minor axis length of
%               crater boundary.
%           Crater_Ellipse: Best-fitting ellipse of crater boundary.
%           Crater_Major_Axis_Azimuth: Azimuthal direction of best-fitting 
%               crater boundary ellipse.
%           Crater_Depth: Depths of crater, defined as distance between 
%               crater minima and surface over minima, given as a structure 
%               which follows the interpSurfaces input.
%           Crater_Maximum_Depth: Maximum depth of crater, defined as 
%               distance between the crater minima and highest boundary 
%               elevation.
%           Crater_Volume: Volumes of crater between crater minima and 
%               crater surface, given as a structure which follows the 
%               interpSurfaces input.
%           Crater_Maximum_Volume: Maximum volume of crater between minima 
%               and highest boundary elevation.
%           Crater_MaxDiam_Ellipticity_Indexes: Ellipticity index of the 
%               crater using the maximum diameter.
%           Crater_BFEllipse_Ellipticity_Indexes: Ellipticity index  
%               of the crater using the best-fitting ellipse major axis.
%           Crater_MaxDiam_Irregularity_Indexes: Irregularity index 
%               of the crater using the maximum diameter.
%           Crater_BFEllipse_Irregularity_Indexes: Irregularity index 
%               of the crater using the best-fitting ellipse major axis.
%           Crater_Mean_Slope: Mean slope within crater.
%           Crater_Median_Slope: Median slope within crater.
%           Crater_Contour_Stats: Contour statistics within the crater,
%               containing investigated contour value, mean slope, median
%                slope, min slope, and max slope.
%           CraterDepth_CraterWidth: Ratio between crater depth and crater
%               width.
%           CraterWidth_BasalWidth: Ratio between crater width and basal
%               width.
%           CraterDepth_BasalHeight: Ratio between crater depth and edifice
%               height.
%           Crater_Axis_Ellipticity: Ellipse axis ellipticity value.
%           Crater_Surface_XY: X-Y coordinates of the crater's surface.
%           Crater_Surface_Z: Surface interpolated elevations for the
%               X-Y coordinates in Crater_Surface_XY, given as a structure 
%               which follows the interpSurfaces input.
%       PeakParams: 
%           Summit_Contour_Peak_Count: Number of contour-based 'peaks' in   
%               the edifice summit region.
%           Main_Flank_Contour_Peak_Count: Number of contour-based 'peaks'   
%               in the edifice main flank region.
%           Lower_Flank_Contour_Peak_Count: Number of contour-based 'peaks'    
%               in the edifice lower flank region.
%           Summit_Local_Peak_Count: Number of local maxima 'peaks' in the 
%               edifice summit region.
%           Main_Flank_Local_Peak_Count: Number of local maxima 'peaks' in   
%               the edifice main flank region.
%           Lower_Flank_Local_Peak_Count: Number of local maxima 'peaks' in    
%               the edifice lower flank region.

%% Unwrap Package
pack = Morvolc_DefaultVals(pack);
disp(sprintf('USING MORVOLC VERSION %s',pack.version))
disp('Unpacking Input and Importing Data...')

morRes = [];
morRes.GeneralParams.inputs = pack;

tifFile = pack.tifFile;
boundaryXYZ = pack.boundaryXY;
contIter = pack.contIter;
maskRegions = pack.maskMap;
dx = pack.dx;
peakDiff = pack.peakDiff;
craterXYZ = pack.craterXY;
craterContIter = pack.craterContIter;
peakContIter = pack.peakContIter;
plotResults = pack.plotResults;
saveFigFolder = pack.saveFigFolder;
saveResFolder = pack.saveResFolder;
figPrefix = pack.figPrefix;
interpSurfaces = pack.interpSurfaces;
xlsFile = pack.xlsFile;
correctCrater = pack.correctCrater;
summitRegion = pack.summitRegion;
ignoreSummitIndex = pack.ignoreSummitIndex;

morRes.GeneralParams.Version = pack.version;
morRes.GeneralParams.StartTime = datetime('now');

Morvolc_CheckIters(contIter,craterXYZ,craterContIter,peakContIter)

%% Setup
% If boundary is given as shapefile, convert to an array.
if ischar(boundaryXYZ)
    Sh = shaperead(boundaryXYZ);
    tx = Sh.X;
    ty = Sh.Y;
    if isnan(tx(end))
        tx = tx(1:end-1);
        ty = ty(1:end-1);
    end
    
    try
        [ttx,tty,~] = ll2utm(ty,tx);
        boundaryXYZ = [ttx',tty'];
    catch
        warning('Unable to convert boundary from Lat/Lon, assuming already in UTM')
        boundaryXYZ = [tx',ty'];
    end
    
    
    boundaryXYZ = double(boundaryXYZ);
end

% If crater is given as shapefile, convert to cell array.
if ischar(craterXYZ)
    Sh = shaperead(craterXYZ);
    craterXYZ = {};
    for i = 1:size(Sh,1)
        tx = Sh(i).X;
        ty = Sh(i).Y;
        if isnan(tx(end))
            tx = tx(1:end-1);
            ty = ty(1:end-1);
        end

        try
            [ttx,tty,~] = ll2utm(ty,tx);
            txyz = [ttx',tty'];
        catch
            warning('Unable to convert crater from Lat/Lon, assuming already in UTM')
            txyz = [tx',ty'];
        end
        
        craterXYZ = [craterXYZ;{double(txyz)}];
    end
elseif ~isempty(craterXYZ) && ~iscell(craterXYZ)
    craterXYZ = {craterXYZ};
end

% If mask given as shapefile, convert to array
maskXY = {};
if ischar(maskRegions)
    Sh = shaperead(maskRegions);
    for i = 1:length(Sh)
        tx = Sh(i).X;
        ty = Sh(i).Y;
        if isnan(tx(end))
            tx = tx(1:end-1);
            ty = ty(1:end-1);
        end

        try
            [ttx,tty,~] = ll2utm(ty,tx);
            maskXY = [maskXY,{[ttx',tty']}];
        catch
            warning('Unable to convert mask region from Lat/Lon, assuming already in UTM')
            maskXY = [maskXY,{double([tx',ty'])}];
        end
    end
elseif isempty(maskRegions)
    maskXY = {};
elseif ~iscell(maskRegions)
    maskXY = {maskRegions};
else
    maskXY = maskRegions;
end

%% Load Tif file, isolate volcano
disp('Loading DEM...')
if ischar(tifFile)
    warning('off','all')
    DEM = GRIDobj(tifFile);
    warning('on','all')
    try
        DEM = reproject2utm(DEM,dx);
    catch
        warning('Cannot reproject DEM to UTM, assuming map is already in coordinates.')
        DEM.Z(DEM.Z<-30000) = NaN;
        DEM = resample(DEM,dx);
    end
else
    DEM = tifFile;
end

morRes.GeneralParams.DEM = DEM;
[Z0,x,y] = GRIDobj2mat(DEM);
[X,Y] = meshgrid(x,y);
X = double(X);
Y = double(Y);

% Crop DEM to boundary 
disp('Clipping Topography...')
bz = interp2(X,Y,Z0,boundaryXYZ(:,1),boundaryXYZ(:,2));
boundaryXYZ = [boundaryXYZ,double(bz)];
boundaryXYZ(sum(isnan(boundaryXYZ),2)>0,:) = [];
p = inpolygon(X,Y,boundaryXYZ(:,1),boundaryXYZ(:,2));
DEM.Z(~p) = NaN;
DEM = crop(DEM);

% Generate X,Y,Z,slope, curvature grids
[Z0,x,y] = GRIDobj2mat(DEM);
Z = Z0;
[X,Y] = meshgrid(x,y);
X = double(X);
Y = double(Y);
Z = double(Z);
S = gradient8(DEM,'deg');
C = curvature(DEM);
[S,~,~] = GRIDobj2mat(S);
[C,~,~] = GRIDobj2mat(C);

% Set everything outside of edifice to NaN.
p = inpolygon(X,Y,boundaryXYZ(:,1),boundaryXYZ(:,2));
Z(~p) = NaN;
S(~p) = NaN;
C(~p) = NaN;

morRes.GeneralParams.X = X;
morRes.GeneralParams.Y = Y;
morRes.GeneralParams.Z = Z;
morRes.GeneralParams.S = S;
morRes.GeneralParams.C = C;

% Correct craters and create crater-masked slope grid
S_nocrater = S;
for i = 1:length(craterXYZ) 
    if correctCrater
        tDEM = fillsinks(DEM,5);
        tF = FLOWobj(tDEM,'preprocess','none');
        dd = drainagebasins(tF);
        [dd,~,~] = GRIDobj2mat(dd);

        mx = mean(craterXYZ{i}(:,1));
        my = mean(craterXYZ{i}(:,2));
        tj = find(abs(x-mx)==min(abs(x-mx)));
        ti = find(abs(y-my)==min(abs(y-my)));
        bb = bwboundaries(dd==dd(ti,tj));
        txyz = [];
        for j = 1:size(bb{1}(:,1))
            txyz = [txyz;x(bb{1}(j,2)),y(bb{1}(j,1))];
        end
    else
        txyz = craterXYZ{i}(:,1:2);
    end

    tz = interp2(X,Y,Z0,txyz(:,1),txyz(:,2));
    craterXYZ{i} = double([txyz,tz]);

    p = inpolygon(X,Y,txyz(:,1),txyz(:,2));
    S_nocrater(p==1) = NaN;
end

%% Determine Contour Interval & Lower Flank Boundary
disp('Determing contours and lower flank boundary...')

if length(contIter) == 1
    if contIter < 0 && contIter > -1
%         trueCont = abs(contIter)*range(Z(:));
        trueCont = round(abs(contIter)*range(Z(:)),0);
    else
%         trueCont = contIter;
        trueCont = round(contIter,0);
    end

    % Determine contours
    tmpConts = round(min(boundaryXYZ(:,3))):round(min(boundaryXYZ(:,3)))+trueCont;
    tmp = mod(tmpConts,trueCont);
    tmpI = find(tmp==0,1);
    conts = [tmpConts(tmpI)-trueCont:trueCont:max(Z(:))+trueCont]';
else
    conts = zeros(length(contIter),1);
    for i = 1:length(contIter)
        if contIter(i) < 0 && contIter(i) > -1
            conts(i) = round(abs(contIter(i))*range(Z(:)),0);
        else
            conts(i) = contIter(i);
        end
    end
    trueCont = mean(diff(conts));
end

% Find Lowest Closed Contour
LowFlankZ = NaN;
LowFlankXYZ = [];
for i = 1:length(conts)
    tmpZ = Z;
    cc = contourc(X(1,:),Y(:,1),Z,[conts(i),conts(i)]);
    cc = Convert_Contours(cc,0);
    tmpZ(tmpZ<conts(i)) = NaN;
    checkArea = sum(~isnan(tmpZ(:)))*dx^2;
    
    for j = 1:length(cc)
        ddx = cc{j}(1,1)-cc{j}(end,1);
        ddy = cc{j}(1,2)-cc{j}(end,2);
        
        if sqrt(ddx^2+ddy^2) <= sqrt(2*(2*dx)^2)
            ar = polyarea(cc{j}(:,1),cc{j}(:,2));
            if ar/checkArea > .5
                LowFlankZ = conts(i);
                LowFlankXYZ = cc{j};
                zz = interp2(X,Y,Z,LowFlankXYZ(:,1),LowFlankXYZ(:,2));
                LowFlankXYZ = [LowFlankXYZ,zz];
                break
            end
        end
    end
    
    if ~isnan(LowFlankZ)
        break
    end
end

if isnan(LowFlankZ)
    LowFlankXYZ = boundaryXYZ;
    LowFlankZ = min(boundaryXYZ(:,3));
end

%% Determine Summit Area
disp('Determing summit area and contour parameters...')

if isempty(craterXYZ)
    useMaxZ = max(Z(:)); 
    if summitRegion==4
        summitRegion = 1;
    end
else
    useMaxZ = Inf;
    for i = 1:length(craterXYZ)
        useMaxZ = min([useMaxZ,min(craterXYZ{i}(:,3))]);
    end
end

[summitXYZ,summitCont] = Morvolc_GetSummitRegion(summitRegion,Z,S,X,Y,LowFlankZ,conts,useMaxZ,peakDiff);

%% Mask grids and boundaries
disp('Applying Mask...')
if ~isempty(maskXY)
    for j = 1:length(maskXY)
        p = inpolygon(X,Y,maskXY{j}(:,1),maskXY{j}(:,2));
        tmpZ = Z;
        Z(p) = NaN;
        S(p) = NaN;
        S_nocrater(p) = NaN;
    
        p = inpolygon(summitXYZ(:,1),summitXYZ(:,2),maskXY{j}(:,1),maskXY{j}(:,2)); 
        if ~p(1)
            tt = find(p == 1,1);
            p = circshift(p,tt-1);
            summitXYZ = circshift(summitXYZ,tt-1,1);
        end
        summitXYZ(p,:) = [];
        
        p = inpolygon(boundaryXYZ(:,1),boundaryXYZ(:,2),maskXY{j}(:,1),maskXY{j}(:,2));
        if ~p(1)
            tt = find(p == 1,1);
            p = circshift(p,tt-1);
            boundaryXYZ = circshift(boundaryXYZ,tt-1,1);
        end
        boundaryXYZ(p,:) = [];
        
        p = inpolygon(LowFlankXYZ(:,1),LowFlankXYZ(:,2),maskXY{j}(:,1),maskXY{j}(:,2));
        if ~p(1)
            tt = find(p == 1,1);
            p = circshift(p,tt-1);
            LowFlankXYZ = circshift(LowFlankXYZ,tt-1,1);
        end
        LowFlankXYZ(p,:) = [];
        
        for i = 1:size(craterXYZ,1)
            tmpC = craterXYZ{i};
            tmpM = maskXY{j};
            p1 = inpolygon(tmpC(:,1),tmpC(:,2),tmpM(:,1),tmpM(:,2));
            p2 = inpolygon(tmpM(:,1),tmpM(:,2),tmpC(:,1),tmpC(:,2));
            [ix,iy] = polyxpoly(tmpC(:,1),tmpC(:,2),tmpM(:,1),tmpM(:,2));
            
            if sum(p1) == 0
                continue;
            end
            
            if sum(p1) == length(p1)
                craterXYZ{i}(p1,:) = [];
                continue;
            end
            
            if p1(1)
                tt = find(p1==0,1);
                p1 = circshift(p1,tt-1);
                tmpC = circshift(tmpC,tt-1,1);
            end
            
            if p2(1)
                tt = find(p2==0,1);
                p2 = circshift(p2,tt-1);
                tmpM = circshift(tmpM,tt-1,1);
            end
            
            p1_t1 = find(p1==1,1);
            p1_t2 = find(p1==1,1,'last');
            p2_t1 = find(p2==1,1);
            p2_t2 = find(p2==1,1,'last');
            
            maskBounds = tmpM(p2_t1:p2_t2,:);
            
            newBound = tmpC(1:p1_t1,1:2);
            dists = sqrt((ix-newBound(end,1)).^2 + (iy-newBound(end,2)).^2);
            useI = find(dists == min(dists));
            newBound = [newBound;ix(useI),iy(useI)];
            
            dists = sqrt((ix(useI)-maskBounds(:,1)).^2 + (iy(useI)-maskBounds(:,2)).^2);
            if dists(end) < dists(1)
                maskBounds = flipud(maskBounds);
            end
            
            newBound = [newBound;maskBounds];
            dists = sqrt((ix-newBound(end,1)).^2 + (iy-newBound(end,2)).^2);
            useI = find(dists == min(dists));
            newBound = [newBound;ix(useI),iy(useI)];
            
            if p1_t2 < size(tmpC,1)
                newBound = [newBound;tmpC(p1_t2+1:end,1:2)];
            end
            
            newZ = interp2(X,Y,tmpZ,newBound(:,1),newBound(:,2));
            craterXYZ{i} = [newBound,newZ];
        end
    end
end

morRes.GeneralParams.MaskXY = maskXY;
morRes.GeneralParams.Lower_Flank_Contour = LowFlankZ;
morRes.GeneralParams.Lower_Flank_XY = LowFlankXYZ;

morRes.GeneralParams.Summit_Contour = summitCont;
morRes.GeneralParams.SummitXYZ = summitXYZ;
morRes.GeneralParams.BoundaryXYZ = boundaryXYZ;
morRes.GeneralParams.CraterXYZ = craterXYZ;

%% Calculate Contour Parameters
disp('Calculating Contour Parameters')

contEllipseAz_ei_ii_meMedSl_minMaxSl_aE = zeros(length(conts),11)*NaN;
contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,1) = conts;
contEls = cell(length(conts),1);
for i = 2:length(conts)-1
    ccXY = contourc(X(1,:),Y(:,1),Z,[conts(i),conts(i)]);
    ccXY = Convert_Contours(ccXY,0);

    if isempty(ccXY)
        continue;
    end
    
    % Isolate Contour to Use for Ellipse
    useJ = 1;
    for j = 2:length(ccXY)
        if size(ccXY{j},1) > size(ccXY{useJ},1)
            useJ = j;
        end
    end
   
    if conts(i) >= LowFlankZ
        ce = FitEllipse_Unwrap(ccXY{useJ});
        if ~isreal(ce.shortAxis)
            contEls{i} = NaN;
            contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,2) = NaN;
            ce.shortAxis = NaN;
            ce.longAxis = NaN;
        else
            contEls{i} = ce;
            contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,2) = ce.phi;
        end
    end
    
    % Contour ellipticity & irregularity indexes.
    if isempty(ccXY) || (ignoreSummitIndex==0 && conts(i) > min(summitXYZ(:,3))) || conts(i) < LowFlankZ
        ei_MaxDiam = NaN;
        ii_MaxDiam = NaN;
        ee = NaN;
        ei_BFEllipse = NaN;
        ii_BFEllipse = NaN;
    else
        [ei_MaxDiam,ii_MaxDiam,~,~] = Morvolc_Ellipticity_Irregularity_MaxDiam(ccXY{useJ}(:,1),ccXY{useJ}(:,2));
        [ei_BFEllipse,ii_BFEllipse,~,~] = Morvolc_Ellipticity_Irregularity_BFEllipse(ccXY{useJ}(:,1),ccXY{useJ}(:,2),ce);
        ee = ce.shortAxis/ce.longAxis;
    end

    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,3) = ei_MaxDiam;
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,4) = ii_MaxDiam;
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,9) = ee;
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,10) = ei_BFEllipse;
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,11) = ii_BFEllipse;
    
    % Get slopes of contour
    ss = S_nocrater;
    ss(Z>conts(i)) = NaN;
    ss(Z<=conts(i-1)) = NaN;
    
    % Fill mean and median slopes.
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,5) = nanmean(ss(:));
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,6) = nanmedian(ss(:));
    
    % Fill min and max Slopes
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,7) = nanmin(ss(:));
    contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(i,8) = nanmax(ss(:));
end    

%% Collect Size Parameters
disp('Collecting Size Parameters...')
% Basal Area & Width
[basalArea,basalWidth] = Morvolc_Calculate_Area_Width(boundaryXYZ(:,1),boundaryXYZ(:,2));

morRes.SizeParams.Basal_Area = basalArea;
morRes.SizeParams.Basal_Width = basalWidth;

% Basal Axes
be = FitEllipse_Unwrap(boundaryXYZ(:,1:2));
morRes.SizeParams.Major_Basal_Axis = be.longAxis;
morRes.SizeParams.Minor_Basal_Axis = be.shortAxis;
morRes.SizeParams.Basal_Axis_Ellipticity = be.shortAxis/be.longAxis;
morRes.SizeParams.Basal_Ellipse = be;

% Summit Area & Width
[summitArea,summitWidth] = Morvolc_Calculate_Area_Width(summitXYZ(:,1),summitXYZ(:,2));

morRes.SizeParams.Summit_Area = summitArea;
morRes.SizeParams.Summit_Width = summitWidth;

% Summit Axes
se = FitEllipse_Unwrap(summitXYZ(:,1:2));
morRes.SizeParams.Summit_Basal_Axis = se.longAxis;
morRes.SizeParams.Summit_Basal_Axis = se.shortAxis;
morRes.SizeParams.Summit_Basal_Axis_Ellipticity = se.shortAxis/se.longAxis;
morRes.SizeParams.Summit_Ellipse = se;

% Heights and Volumes
[surfXYs,surfZs,maxHeight,maxVol,volcHeights,volcVols,anyHeights] = Morvolc_CalculateHeightVols(X,Y,Z,boundaryXYZ,min(boundaryXYZ(:,3)),interpSurfaces);

morRes.SizeParams.Peak_Height = volcHeights;
morRes.SizeParams.Maximum_Height = maxHeight;
morRes.SizeParams.Any_Height = anyHeights;
morRes.GeneralParams.Basal_Surface_XY = surfXYs;
morRes.GeneralParams.Basal_Surface_Z = surfZs;

morRes.SizeParams.Volume = volcVols;
morRes.SizeParams.Maximum_Volume = maxVol;

% Eroded Volume
[cont_area_convHullArea,volumeDiff,ConvZ] = CalculateErodedVolume(X,Y,Z,trueCont,min(boundaryXYZ(:,3)),1,1);
morRes.SizeParams.Minimum_Eroded_Volume = volumeDiff;
morRes.SizeParams.Convex_Hull_Areas = cont_area_convHullArea;
morRes.SizeParams.Convex_Hull_Interpolated_Surface = ConvZ;

%% Collect Orientation Parameters
disp('Collecting Orientation Parameters...')
morRes.OrientParams.Basal_Major_Axis_Azimuth = be.phi;
morRes.OrientParams.Summit_Major_Axis_Azimuth = se.phi;
morRes.OrientParams.Contour_Values = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,1);
morRes.OrientParams.Contour_Major_Axis_Azimuths = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,2);

%% Collect Shape Parameters
disp('Collecting Shape Parameters...')
morRes.ShapeParams.Contour_Values = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,1);
morRes.ShapeParams.Contour_Axis_Ellipticity = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,9);

morRes.ShapeParams.Contour_MaxDiam_Ellipticity_Indexes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,3);
morRes.ShapeParams.Contour_MaxDiam_Irregularity_Indexes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,4);
morRes.ShapeParams.Mean_MaxDiam_Ellipticity_Index = nanmean(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,3));
morRes.ShapeParams.Mean_MaxDiam_Irregularity_Index = nanmean(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,4));


morRes.ShapeParams.Contour_BFEllipse_Ellipticity_Indexes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,10);
morRes.ShapeParams.Contour_BFEllipse_Irregularity_Indexes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,11);
morRes.ShapeParams.Mean_BFEllipse_Ellipticity_Index = nanmean(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,10));
morRes.ShapeParams.Mean_BFEllipse_Irregularity_Index = nanmean(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,11));

morRes.ShapeParams.Mean_Ellipse_Ellipticity = nanmean(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,9));
morRes.ShapeParams.Height_BasalWidth.Natural = morRes.SizeParams.Peak_Height.Natural/morRes.SizeParams.Basal_Width;
morRes.ShapeParams.Height_BasalWidth.IDW = morRes.SizeParams.Peak_Height.IDW/morRes.SizeParams.Basal_Width;
morRes.ShapeParams.Height_BasalWidth.Kringing = morRes.SizeParams.Peak_Height.Kringing/morRes.SizeParams.Basal_Width;
morRes.ShapeParams.SummitWidth_BasalWidth = morRes.SizeParams.Summit_Width/morRes.SizeParams.Basal_Width;
morRes.ShapeParams.Contour_Ellipses = contEls;

Znan = isnan(Z);
Z(isnan(Z)) = 0;
[~,~,~,~,~,~,sk,ku] = Calculate_Moments(x,y',Z);
Z(Znan) = NaN;
morRes.ShapeParams.Skewness = sk;
morRes.ShapeParams.Kurtosis = ku;

%% Collect Slope Parameters
disp('Collecting Slope Parameters...')
tmpS = S_nocrater;
tmpS(isnan(Z)) = NaN;
allAvgSl = nanmean(tmpS(:));
allMedSl = nanmedian(tmpS(:));
allStdSl = nanstd(tmpS(:));
morRes.SlopeParams.WholeEdifice_Mean_Slope = allAvgSl;
morRes.SlopeParams.WholeEdifice_Median_Slope = allMedSl;
morRes.SlopeParams.WholeEdifice_Std_Slope = allStdSl;

pF = inpolygon(X,Y,LowFlankXYZ(:,1),LowFlankXYZ(:,2));
pS = inpolygon(X,Y,summitXYZ(:,1),summitXYZ(:,2));

tmpS = S_nocrater;
tmpS(pF) = NaN;
lFlankAvgSl = nanmean(tmpS(:));
lFlankMedSl = nanmedian(tmpS(:));
lFlankStdSl = nanstd(tmpS(:));
morRes.SlopeParams.Lower_Flank_Mean_Slope = lFlankAvgSl;
morRes.SlopeParams.Lower_Flank_Median_Slope = lFlankMedSl;
morRes.SlopeParams.Lower_Flank_Std_Slope = lFlankStdSl;

tmpS = S_nocrater;
tmpS(~pS) = NaN;
summitAvgSl = nanmean(tmpS(:));
summitMedSl = nanmedian(tmpS(:));
summitStdSl = nanstd(tmpS(:));
morRes.SlopeParams.Summit_Mean_Slope = summitAvgSl;
morRes.SlopeParams.Summit_Median_Slope = summitMedSl;
morRes.SlopeParams.Summit_Std_Slope = summitStdSl;

tmpS = S_nocrater;
tmpS(~pF) = NaN;
tmpS(pS) = NaN;
flankAvgSl = nanmean(tmpS(:));
flankMedSl = nanmedian(tmpS(:));
flankStdSl = nanstd(tmpS(:));
morRes.SlopeParams.Flank_Mean_Slope = flankAvgSl;
morRes.SlopeParams.Flank_Median_Slope = flankMedSl;
morRes.SlopeParams.Flank_Std_Slope = flankStdSl;

morRes.SlopeParams.Contour_Values = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,1);
morRes.SlopeParams.Contour_Mean_Slopes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,5);
morRes.SlopeParams.Contour_Median_Slopes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,6);
morRes.SlopeParams.Contour_Min_Slopes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,7);
morRes.SlopeParams.Contour_Max_Slopes = contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,8);
morRes.SlopeParams.Contour_Max_Mean_Slope = max(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,5));
morRes.SlopeParams.Contour_Max_Median_Slope = max(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,6));

hf_i = find(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,5)==max(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,5)));
hf_cont= contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(hf_i,1);
totHeight = max(Z(:)) - max(boundaryXYZ(:,3));
heightFrac = (hf_cont-max(boundaryXYZ(:,3)))/totHeight;

morRes.SlopeParams.Max_Slope_Height = hf_cont;
morRes.SlopeParams.Max_Slope_Height_Fraction = heightFrac;
morRes.SlopeParams.Max_Slope = max(contEllipseAz_ei_ii_meMedSl_minMaxSl_aE(:,5));

%% Collect Crater Parameters
disp('Collecting Crater Parameters...')
if isempty(craterXYZ)
    morRes.CraterParams = [];
    morRes.GeneralParams.Crater_Interp = [];
else
    craterArea = zeros(size(craterXYZ))*NaN;
    craterWidth = craterArea;
    craterMajor = craterArea;
    craterMinor = craterArea;
    craterEllipse = cell(size(craterXYZ));
    craterAz = craterArea;
    craterDepth = craterEllipse;
    craterMaxDepth = craterArea;
    craterVol = craterEllipse;
    craterSurfXY = craterEllipse;
    craterSurfZ = craterEllipse;
    craterMaxVol = craterArea;
    craterElipt_MaxDiam = craterArea;
    craterIreg_MaxDiam = craterArea;
    craterElipt_BFEllipse = craterArea;
    craterIreg_BFEllipse = craterArea;
    craterMeanSl = craterArea;
    craterMedianSl = craterArea;
    craterContStats = craterEllipse;
    craterDepth_Width = craterEllipse;
    craterWidth_BasalWidth = craterArea;
    craterDepth_BasalHeight = craterEllipse;
    craterAxisEllipticity = craterArea;
    
    for i = 1:length(craterXYZ)
        cxyz = craterXYZ{i};
        
        % Area & Width
        [ca,cw] = Morvolc_Calculate_Area_Width(cxyz(:,1),cxyz(:,2));
        
        craterArea(i) = ca;
        craterWidth(i) = cw;
        
        % Ellipse
        cEl = FitEllipse_Unwrap(cxyz(:,1:2));
        
        craterMajor(i) = cEl.longAxis;
        craterMinor(i) = cEl.shortAxis;
        craterAz(i) = cEl.phi;
        craterEllipse{i} = cEl;
        craterAxisEllipticity(i) = cEl.shortAxis/cEl.longAxis;
        
        % Irregularity & Ellipticity
        [cEi_MaxDiam,cIi_MaxDiam,~,~] = Morvolc_Ellipticity_Irregularity_MaxDiam(cxyz(:,1),cxyz(:,2));
        [cEi_BFEllipse,cIi_BFEllipse,~,~] = Morvolc_Ellipticity_Irregularity_BFEllipse(cxyz(:,1),cxyz(:,2),cEl);
        
        craterElipt_MaxDiam(i) = cEi_MaxDiam;
        craterIreg_MaxDiam(i) = cIi_MaxDiam;
        craterElipt_BFEllipse(i) = cEi_BFEllipse;
        craterIreg_BFEllipse(i) = cIi_BFEllipse;
        
        % Depth & Volume
        tmpZ = Z;
        p = inpolygon(X,Y,cxyz(:,1),cxyz(:,2));
        tmpZ(~p) = NaN;
        tmpCxyz = cxyz;
        [cSurfXYs,cSurfZs,mH,mV,vHs,vVs] = Morvolc_CalculateHeightVols_Crater(X,Y,tmpZ,tmpCxyz,max(tmpCxyz(:,3)),interpSurfaces);

        craterMaxDepth(i) = mH;
        craterMaxVol(i) = mV;
        craterDepth{i} = vHs;
        craterVol{i} = vVs;
        craterSurfXY{i} = cSurfXYs;
        craterSurfZ{i} = cSurfZs;
  
        % Slopes
        tmpS = S;
        tmpS(~p) = NaN;
        craterMeanSl(i) = nanmean(tmpS(:));
        craterMedianSl(i) = nanmedian(tmpS(:));

        tmpZ = Z;
        tmpZ(~p) = NaN;
        
        % Contour Stats
        if craterContIter < 0 && craterContIter > -1
%             trueCraterCont = round(abs(craterContIter)*range(tmpZ(:)),0);
            trueCraterCont = round(abs(craterContIter)*range(tmpZ(:)));
        else
%             trueCraterCont = round(craterContIter,0);
            trueCraterCont = round(craterContIter);
        end
        
        tmpConts = round(min(tmpZ(:))):round(min(tmpZ(:)))+trueCraterCont;
        tmp = mod(tmpConts,trueCraterCont);
        tmpI = find(tmp==0,1);
        craterConts = [tmpConts(tmpI)-trueCraterCont:trueCraterCont:max(tmpZ(:))]';
        
        
        tmpStats.Crater_Contours = craterConts;
        tmpCMeanSl = zeros(size(tmpStats.Crater_Contours));
        tmpCMedianSl = zeros(size(tmpStats.Crater_Contours));
        tmpCMinSl = zeros(size(tmpStats.Crater_Contours));
        tmpCMaxSl = zeros(size(tmpStats.Crater_Contours));
        
        for j = 2:length(craterConts)
            ccI = contourc(tmpZ,[craterConts(j),craterConts(j)]);
            ccI = Convert_Contours(ccI,0);

            ss = tmpS;
            ss(tmpZ>craterConts(j)) = NaN;
            ss(tmpZ<=craterConts(j-1)) = NaN;
            
            tmpCMeanSl(j) = nanmean(ss(:));
            tmpCMedianSl(j) = nanmedian(ss(:));
            tmpCMinSl(j) = nanmin(ss(:));
            tmpCMaxSl(j) = nanmax(ss(:));
        end
        
        tmpStats.Crater_Mean_Slope = tmpCMeanSl;
        tmpStats.Crater_Median_Slope = tmpCMedianSl;
        tmpStats.Crater_Min_Slope = tmpCMinSl;
        tmpStats.Crater_Max_Slope = tmpCMaxSl;
        
        craterContStats{i} = tmpStats;
        
        % Other Crater Stats
        craterDepth_Width{i}.Natural = craterDepth{i}.Natural/craterWidth(i);
        craterDepth_Width{i}.IDW = craterDepth{i}.IDW/craterWidth(i);
        craterDepth_Width{i}.Kringing = craterDepth{i}.Kringing/craterWidth(i);
        craterWidth_BasalWidth(i) = craterWidth(i)/basalWidth;
        craterDepth_BasalHeight{i}.Natural = craterDepth{i}.Natural/volcHeights.Natural;
        craterDepth_BasalHeight{i}.IDW = craterDepth{i}.IDW/volcHeights.IDW;
        craterDepth_BasalHeight{i}.Kringing = craterDepth{i}.Kringing/volcHeights.Kringing;
    end
    
    morRes.CraterParams.Crater_Area = craterArea;
    morRes.CraterParams.Crater_Width = craterWidth;
    morRes.CraterParams.Crater_Major_Axis = craterMajor;
    morRes.CraterParams.Crater_Minor_Axis = craterMinor;
    morRes.CraterParams.Crater_Ellipse = craterEllipse;
    morRes.CraterParams.Crater_Major_Axis_Azimuth = craterAz;
    morRes.CraterParams.Crater_Depth = craterDepth;
    morRes.CraterParams.Crater_Maximum_Depth = craterMaxDepth;
    morRes.CraterParams.Crater_Volume = craterVol;
    morRes.CraterParams.Crater_Maximum_Volume = craterMaxVol;
    morRes.CraterParams.Crater_BFEllipse_Ellipticity_Index = craterElipt_BFEllipse;
    morRes.CraterParams.Crater_BFEllipse_Irregularity_Index = craterIreg_BFEllipse;
    morRes.CraterParams.Crater_MaxDiam_Ellipticity_Index = craterElipt_MaxDiam;
    morRes.CraterParams.Crater_MaxDiam_Irregularity_Index = craterIreg_MaxDiam;
    morRes.CraterParams.Crater_Mean_Slope = craterMeanSl;
    morRes.CraterParams.Crater_Median_Slope = craterMedianSl;
    morRes.CraterParams.Crater_Contour_Stats = craterContStats;
    morRes.CraterParams.CraterDepth_CraterWidth = craterDepth_Width;
    morRes.CraterParams.CraterWidth_BasalWidth = craterWidth_BasalWidth;
    morRes.CraterParams.CraterDepth_BasalHeight = craterDepth_BasalHeight;
    morRes.CraterParams.Crater_Axis_Ellipticity = craterAxisEllipticity;
    morRes.CraterParams.Crater_Surface_XY = craterSurfXY;
    morRes.CraterParams.Crater_Surface_Z = craterSurfZ;
end

%% Collect Peak Parameters
disp('Collecting Peak Parameters...')

if peakContIter < 0 && peakContIter > -1
%     truePeakCont = round(abs(peakContIter)*maxHeight,0);
    truePeakCont = abs(peakContIter)*maxHeight;
else
%     truePeakCont = round(peakContIter,0);
    truePeakCont = peakContIter;
end

allPeakIJ_cont = [];
conts = min(Z(:)):truePeakCont:max(Z(:));

% Contour-based Peak Count
for i = length(conts):-1:1
    bb = bwboundaries(Z>=conts(i));
    for j = 1:length(bb)
        if ~isempty(allPeakIJ_cont)
            [in,on] = inpolygon(allPeakIJ_cont(:,2),allPeakIJ_cont(:,1),bb{j}(:,2),bb{j}(:,1));
            ion = in+on;

            if sum(ion) == 0
                allPeakIJ_cont = [allPeakIJ_cont;bb{j}(1,:)];
            end
        else
            allPeakIJ_cont = [allPeakIJ_cont;bb{j}(1,:)];
        end
    end
end

allPeakXY_cont = zeros(size(allPeakIJ_cont));
for i = 1:size(allPeakIJ_cont,1)
    allPeakXY_cont(i,:) = [x(allPeakIJ_cont(i,2)),y(allPeakIJ_cont(i,1))];
end

[sinp,son] = inpolygon(allPeakXY_cont(:,1),allPeakXY_cont(:,2),summitXYZ(:,1),summitXYZ(:,2));
[lin,lon] = inpolygon(allPeakXY_cont(:,1),allPeakXY_cont(:,2),LowFlankXYZ(:,1),LowFlankXYZ(:,2));
morRes.PeakParams.Summit_Contour_Peak_Count = sum((sinp+son)>0);
morRes.PeakParams.Main_Flank_Contour_Peak_Count = sum((lin+lon)>0)-sum((sinp+son)>0);
morRes.PeakParams.Lower_Flank_Contour_Peak_Count = size(allPeakXY_cont,1)-sum((lin+lon)>0);

% Local Maxima Peak Count
locMax1 = islocalmax(Z,1);
locMax2 = islocalmax(Z,2);
locMax3 = locMax1.*locMax2;

lM3_summit = locMax3;
[sinp,son] = inpolygon(X,Y,summitXYZ(:,1),summitXYZ(:,2));
lM3_summit((sinp+son)*1 ==0) = NaN;

lM3_lFlank = locMax3;
[lin,lon] = inpolygon(X,Y,LowFlankXYZ(:,1),LowFlankXYZ(:,2));
lM3_lFlank((lin+lon)*1 ==0) = NaN;

morRes.PeakParams.Summit_Local_Peak_Count = nansum(lM3_summit(:));
morRes.PeakParams.Main_Flank_Local_Peak_Count = nansum(lM3_lFlank(:))-nansum(lM3_summit(:));
morRes.PeakParams.Lower_Flank_Local_Peak_Count = nansum(locMax3(:))-nansum(lM3_lFlank(:));

%% Save Results
morRes.GeneralParams.EndTime = datetime('now');
if ~isempty(saveResFolder)
    disp('Saving Results...')
    save([saveResFolder,figPrefix,'MorVolc_Results.mat'],'morRes')
end

if ~isempty(xlsFile)
    Morvolc_SaveXLS(xlsFile,morRes);
end

%% Plot Results
if plotResults
    disp('Plotting Results...')
    MorVolc_Plots(morRes);
end

end
