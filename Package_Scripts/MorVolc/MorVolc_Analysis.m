function MV_Res = MorVolc_Analysis(pack)
%%
% Name: MorVolc_Analysis
% Author: Daniel O'Hara
% Original MorVolc Author: Pablo Grosse
% Initial Date: 02/24/2021 (mm/dd/yyyy)
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
%       tifFile: .tif file to analyze. This can also be given as a
%           pre-compiled GRIDobj object.
%       boundaryXY: 2xN matrix or shapefile name of edifice boundary 
%           locations. If matrix, first column is X and second column is Y.
%       maskXY: Array, cell array, or .shp file and path, of mask regions  
%           X- and Y-coordinates to ignore in analysis.
%       craterXY: 2xN matrix, cell array, or shapefile name of crater 
%           boundary location that is ignored during DEM filling. If 
%           matrix, first column is X and second column is Y.
%
%       dx: Resolution of DEM (in meters).
%
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
%               4: Summit defined by the lowest elevation of the crater 
%                   outline.
%               5: Summit defined by the lowest elevation of the crater.
%               Value raging -1 to 0: Summit defined by upper percentile of
%                   edifice height (similar to DrainageVolc basinTopN flag.
%       ignoreSummitIndex: Flag to ignore summit and perform ellipticity
%           and irregularity analysis over edifice.
%
%       hypsIter: Iteration value (0-1) to use for hypsometric/CDF 
%           analysis.
%       roughnessWindows: Window sizes (in m) to calculated roughness.
%       roughnessType: String describing the type of roughness analysis,
%           values given in TopoToolbox's 'roughness' script.
%       slopeVarianceWindows: Window sizes (in m) to calculate slope
%           variance.
%
%       interpSurfaces: Structure to calculate the basal and crater
%           surfaces using a variety of interpolation values, including:
%               Natural: Matlab natural interpolation.
%               IDW: Inverse Distance Weighting interpolation.
%               Kriging: Kriging interpolation.
%
%       xlsFile: MS Excel file to save analysis results. Leave empty to not
%           save values.
%       saveFigFolder: Folder to save analysis plots (must include '/' or 
%           '\' at the end, and plotResults must be set to 1 to use). Set 
%           to empty if not saving plots.
%       saveInputs: Flag to save the inputs as a text file.
%
%       plotResults: Flag to plot analysis results.
%       visPlots: Flag to determine whether plots are visible (good for
%           running in background).
%       figPrefix: Name prefix for all figures.
%       figTitlePrefix: Prefix for figure titles (e.g, volcano name).
%       saveResFolder: Folder to save output structure. Set to empty if not
%           saving results.
%
%       verbose: Flag to indicate the amount of output displayed in the
%           terminal. 0 gives no output; 1 gives broad output of each step;
%           2 gives detailed information.
%       zipFiles: Flag to determine whether the folders contain the results
%           and or figures should be automatically zipped (helpful for
%           running scripts on HPCs or for large amounts of data). 1 will 
%           zip only the figures folder, 2 will zip only the results 
%           folder, 3 will zip both. If the figure and result folders are 
%           the same, everything is zipped into one folder.  
%       deleteAfterZip: Flag to determine of folders should be
%           automaticallly deleted after being zipped (helpful for running 
%           scripts on HPC or for large amounts of data). 1 will delete the
%           figures folder after zipping, 2 will delete the results folder
%           after zipping, 3 will delete both after zipping.
%
% Output:
%   morRes: Analysis results given as a structure. Fields are:
%       GeneralParams:
%           inputs: Script input structure.
%           Version: Analysis script version number.
%           StartTime: DateTime object giving the time the script was
%               started.
%           EndTime: DateTime object giving the time the script ended.
%
%       GeographicParams:
%           DEM0: Full GRIDobj of DEM (used by TopoToolbox).
%           DEM: GRIDobj of DEM cut to the edifice.
%
%           boundaryXY: Volcano boundary X- and Y-coordinates.
%           craterXY: Crater region X- and Y-coordinates.
%           maskXY: Mask region x- and y-coordinates.
%
%           Slope: Grid of DEM slope values.
%           Lower_Flank_Contour: Elevation between main flank and lower
%               flank.
%           Lower_Flank_XY: X-Y coordinates of the lower flank boundary.
%           Summit_Contour: Elevation between main flank and summit.
%           SummitXYZ: Matrix of summit boundary x, y, and z values.
%           Basal_Surfaces: Structure of basal surfaces derived using
%               the Natural, IDW, and Kriging methods (based on user
%               input). Given as GRIDobjs.
%
%       TopoParams:
%           Hypsometry: Structure of hypsometry values:
%               Elevation_Hyps_Area_Values: mx2 matrix of elevation-based
%                   hypsometry values. Second row is elevations, first row
%                   is associated area.
%               Slope_Hyps_Area_Values: mx2 matrix of slope-based
%                   hypsometry values. Second row is slope, first row
%                   is associated area.
%               ProfileC_Hyps_Area_Values: mx2 matrix of profile curvature
%                   hypsometry values. Second row is profile curvature, 
%                   first row is associated area.
%               PlanformC_Hyps_Area_Values: mx2 matrix of planform
%                   curvature hypsometry values. Second row is planform 
%                   curvature, first row is associated area.
%           SlopeVariance: Structure of slope variance values:
%               Elevation: Structure of slope variance of
%                       elevational bins along the edifice.
%                   Values: Matrix of values, columns are elevation bin, 
%                       normalized elevation bin, mean slope, std slope,  
%                       and slope variance
%                   Titles: Cell array describing Values matrix.
%               Windows: Indexed-structure of slope variance analysis:
%                   SlopeVariance_Grid = Grid of roughness values.
%                   WindowSize = Pixel-size of roughness window.
%                   ExpectedWindowRes = Expected, user-defined roughness
%                       window.
%                   TrueWindowRes = Actual roughness based on grid 
%                       resolution.
%                   Hypsometry_Areas = Normalized area values for 
%                       roughness.
%                   Hypsometry_Values = Roughness CDF values.
%                   Hypsometry_NormValues = Nomalized roughness CDF values.
%               Total: Total slope variance of the entire edifice.
%           Roughness: Indexed-structure of roughness analysis.
%               Roughness_Grid = Grid of roughness values.
%               WindowSize = Pixel-size of roughness window.
%               ExpectedWindowRes = Expected, user-defined roughness
%                   window.
%               TrueWindowRes = Actual roughness based on grid resolution.
%               Hypsometry_Areas = Normalized area values for roughness.
%               Hypsometry_Values = Roughness CDF values.
%               Hypsometry_NormValues = Nomalized roughness CDF values.
%           Curvature_Grids: Structure of curvature GRIDobjs:
%               Profile: Profile curvature GRIDobj.
%               Planform: Planform curvature GRIDobj.
%           Centers: Structure of center locations.
%               Topographic_XY: X-Y location of highest topography.
%               Geometric_XY: X-Y location of center of boundary.
%               Volumetric_XY: X-Y location of center of volume.
%
%       SizeParams:
%           Basal: Structure of parameters at the edifice base.
%               Area: Area of edifice boundary.
%               Width: Diameter of circle with area equivilant to
%                   edifice boundary.
%               Major_Axis: Best-fitting ellipse major axis length of
%                   edifice boundary.
%               Minor_Axis: Best-fitting ellipse minor axis length of
%                   edifice boundary.
%               Axis_Ellipticity: Ellipticity (minor / major axis) value
%                   of the edifice boundary.
%               Ellipse: Best-fitting ellipse of edifice boundary.
%           Summit: Struvture of parameters at the edifice summit.
%               Area: Area of summit boundary.
%               Width: Diameter of circle with area equivilant to
%                   summit boundary.
%               Major_Axis: Best-fitting ellipse major axis length of
%                   summit boundary.
%               Minor_Axis: Best-fitting ellipse minor axis length of
%                   summit boundary.
%               Axis_Ellipticity: Ellipticity (minor / major axis) value
%                   of the summit region.
%               Ellipse:Best-fitting ellipse of summit boundary.
%           Heights: Structure of edifice heights:
%               Surface_to_Peak: Heights of edifice, defined as distance  
%                   between peak and basal surface under peak, given as a  
%                   structure based on interpSurfaces input.
%               Maximum_From_BasalSurface: Maximum edifice height between 
%                   the edifice and basal surface, given as a structure 
%                   based on interpSurfaces input.
%               Maximum_From_Boundary: Maximum height of edifice from the 
%                   lowest edifice boundary to the peak.
%           Volumes: Structure of edifice volumes:
%               Total: Volumes of edifice between peak and basel surface, 
%                   given as a structure based on interpSurfaces input.
%               Maximum: Maximum volume of edifice between peak and
%                   lowest boundary elevation.
%               Minimum_Eroded_Volume: Minimum edifice eroded volume 
%                   calculated as the integrated areal difference between 
%                   actual topography and topography fit by a convex hull 
%                   over contours (O'Hara & Karlstrom, 2023).
%           ReconstructedTopo: Structure of reconstructed topography from
%                   the convex-hull method of O'Hara & Karlstrom, 2023.
%               Convex_Hull_Areas: Contour and convex hull contour areas 
%                   from the eroded volume analysis.
%               Convex_Hull_Interpolated_Surface: Grid of reconstructed
%                   elevations from the convex hull algorithm.
%           
%       OrientParams:
%           Basal_Major_Axis_Azimuth: Azimuthal direction of best-fitting
%               edifice boundary ellipse.
%           Summit_Major_Axis_Azimuth: Azimuthal direction of best-fitting
%               edifice summit ellipse.
%           Contour_Elevation_Major_Axis_Azimuth: mx2 matrix of contour 
%               elevations (first row) and azimuthal direction of 
%               best-fitting ellipses along each contour (second row).
%
%       ShapeParams:
%           Contour: Structure of contour-based values.
%               Elevations: Contour values used for contour analysis.
%               Ellipses: Best-fitting ellipses along each contour.
%               Axis_Ellipticity: Axis ellipticity indexes
%                   (best-fitting ellipse minor axis / major axis).
%               Ellipticity_Indices: Structure of ellipticity indices 
%                       broken down by methodology:
%                   MaxDiameter: Ellipticity index using the maximum 
%                       diameter of the contour (from Grosse et al., 2012).
%                   BFEllipse: Ellipticity index using the contour's 
%                       best-fitting ellipse major axis.
%               Irregularity_Indices: Structure of irregularity indices 
%                       broken down by methodology:
%                   MaxDiameter: Irregularity index using the maximum 
%                       diameter of the contour (from Grosse et al., 2012).
%                   BFEllipse: Irregularity index using the contour's 
%                       best-fitting ellipse major axis.
%           Mean_Shapes: Structure of average values over all contours:
%               Ellipticity_Index: Structure of mean contour ellipticity 
%                       index broken down by methodology:
%                   MaxDiameter: Mean value of max diameter indices.
%                   BFEllipse: Mean value of best-fitting ellipse indices.
%               Irregularity_Index: Structure of mean contour irregularity 
%                       index broken down by methodology:
%                   MaxDiameter: Mean value of max diameter indices.
%                   BFEllipse: Mean value of best-fitting ellipse indices.
%               Ellipse_Ellipticity: Mean ellipse ellipticity (best-fitting 
%                   ellipse minor axis / major axis).
%           Geometry: Structure of geometric values:
%               Height_BasalWidth: Ratio between edifice height and basal
%                   width, seperated by edifice height derived from basal
%                   surface as given by user input.
%               SummitWidth_BasalWidth: Ratio between summit width and 
%                   basal width.
%           Skewness: Edifice X-Y skewness.
%           Kurtosis: Edifice X-Y kurtosis.
%
%       SlopeParams:
%           Full_Edifice: Structure of slope values considering the entire
%                   edifice:
%               Mean: Mean slope.
%               Median: Median slope.
%               Std: Slope standard deviation.
%           Lower_Flank: Structure of slope values considering the lower
%                   edifice flank (region between lowest and highest 
%                   edifice boundaries):
%               Mean: Mean slope.
%               Median: Median slope.
%               Std: Slope standard deviation.
%           Main_Flank: Structure of slope values considering the main
%                   edifice flank (region between highest edifice boundary 
%                   and summit):
%               Mean: Mean slope.
%               Median: Median slope.
%               Std: Slope standard deviation.
%           Summit: Structure of slope values considering the edifice
%                   summmit:
%               Mean: Mean slope.
%               Median: Median slope.
%               Std: Slope standard deviation.
%           Contour: Structure of contour slope values.
%               Elevations: Contour elevations.
%               Mean: Mean slope along each contour.
%               Median: Median slope along each contour.
%               Min: Min slope along each contour.
%               Max: Max slope along each contour.
%               Std: Slope standard deviation along each contour.
%               Maximum_Mean: Maximum mean slope of all contours.
%               Maximum_Median: Maximum median slope of all contours.
%           Maximum_Slope: Structure of information related to maximum
%                   contour slope.
%               Elevation: Elevation (above sea level) maximum contour  
%                   slope.
%               Height: Height (from edifice base) of the maximum contour 
%                   slope.
%               Height_Fraction: Percentage of maximum contour slope height
%                   compared to entire edifice relief. 
%               Slope: Maximum mean contour slope value.
%
%       RoughnessParams:
%           Full_Edifice: Structure of roughness value parameters for the 
%                   entire edifice.
%               Mean: 1xn array of mean roughness values, where n is the
%                   number of provided roughness windows.
%               Median: 1xn array of median roughness values, where n is 
%                   the number of provided roughness windows.
%               Std: 1xn array of roughness value standard deviations, 
%                   where n is the number of provided roughness windows.
%           Lower_Flank: Structure of roughness value parameters for the 
%                   lower edifice flank (region between the lowest and 
%                   highest edifice boundary).
%               Mean: 1xn array of mean roughness values, where n is the
%                   number of provided roughness windows.
%               Median: 1xn array of median roughness values, where n is 
%                   the number of provided roughness windows.
%               Std: 1xn array of roughness value standard deviations, 
%                   where n is the number of provided roughness windows.
%           Main_Flank: Structure of roughness value parameters for the 
%                   main edifice flank (region between the highest edifice 
%                   boundary and summit).
%               Mean: 1xn array of mean roughness values, where n is the
%                   number of provided roughness windows.
%               Median: 1xn array of median roughness values, where n is 
%                   the number of provided roughness windows.
%               Std: 1xn array of roughness value standard deviations, 
%                   where n is the number of provided roughness windows.
%           Summit: Structure of roughness value parameters for the 
%                   edifice summit.
%               Mean: 1xn array of mean roughness values, where n is the
%                   number of provided roughness windows.
%               Median: 1xn array of median roughness values, where n is 
%                   the number of provided roughness windows.
%               Std: 1xn array of roughness value standard deviations, 
%                   where n is the number of provided roughness windows.
%           Contour: Structure of roughness values for each contour:
%               Elevations: Contour elevations.
%               Mean: Mean roughness along each contour for each window.
%               Median: Median roughness along each contour for each 
%                   window.
%               Min: Min roughness along each contour for each window.
%               Max: Max roughness along each contour for each window.
%               Std: Roughness standard deviation along each contour.
%               Maximum_Mean: Maximum mean roughness of all contours for 
%                   each window.
%               Maximum_Median: Maximum median roughness of all contours 
%                   for each window.
%
%       WindowedSlopeVarianceParams:
%           Full_Edifice: Structure of windowed slope variance (WSV) value 
%                   parameters for the entire edifice.
%               Mean: 1xn array of mean WSV values, where n is the
%                   number of provided WSV windows.
%               Median: 1xn array of median WSV values, where n is 
%                   the number of provided WSV windows.
%               Std: 1xn array of WSV value standard deviations, 
%                   where n is the number of provided WSV windows.
%           Lower_Flank: Structure of WSV value parameters for the 
%                   lower edifice flank (region between the lowest and 
%                   highest edifice boundary).
%               Mean: 1xn array of mean WSV values, where n is the
%                   number of provided WSV windows.
%               Median: 1xn array of median WSV values, where n is 
%                   the number of provided WSV windows.
%               Std: 1xn array of WSV value standard deviations, 
%                   where n is the number of provided WSV windows.
%           Main_Flank: Structure of WSV value parameters for the 
%                   main edifice flank (region between the highest edifice 
%                   boundary and summit).
%               Mean: 1xn array of mean WSV values, where n is the
%                   number of provided WSV windows.
%               Median: 1xn array of median WSV values, where n is 
%                   the number of provided WSV windows.
%               Std: 1xn array of WSV value standard deviations, 
%                   where n is the number of provided WSV windows.
%           Summit: Structure of WSV value parameters for the 
%                   edifice summit.
%               Mean: 1xn array of mean WSV values, where n is the
%                   number of provided WSV windows.
%               Median: 1xn array of median WSV values, where n is 
%                   the number of provided WSV windows.
%               Std: 1xn array of WSV value standard deviations, 
%                   where n is the number of provided WSV windows.
%           Contour: Structure of WSV values for each contour:
%               Elevations: Contour elevations.
%               Mean: Mean WSV along each contour for each window.
%               Median: Median WSV along each contour for each window.
%               Min: Min WSV along each contour for each window.
%               Max: Max WSV along each contour for each window.
%               Std: WSV standard deviation along each contour.
%               Maximum_Mean: Maximum mean WSV of all contours for 
%                   each window.
%               Maximum_Median: Maximum median WSV of all contours 
%                   for each window.
%
%       CraterParams:
%           Size: Structure of crater size parameters, similar to SizeParam
%                   structure:
%               Area: Area of crater boundary.
%               Width: Diameter of circle with area equivilant to
%                   crater boundary.
%               Major_Axis: Best-fitting ellipse major axis length of
%                   crater boundary.
%               Minor_Axis: Best-fitting ellipse minor axis length of
%                   crater boundary.
%               Axis_Ellipticity: Ellipse axis ellipticity value.
%               Ellipse: Best-fitting ellipse of crater boundary.
%               Depth: Structure of crater depth parameters:
%                   Total: Depths of crater, defined as distance between 
%                       crater minima and surface over minima, given as a  
%                       structure which follows the interpSurfaces input.
%                   Maximum: Maximum depth of crater, defined as 
%                       distance between the crater minima and highest  
%                       boundary elevation.
%               Volume: Structure of crater depth parameters:
%                   Total: Volumes of crater between crater minima and 
%                       crater surface, given as a structure which follows
%                       the interpSurfaces input.
%                   Maximum: Maximum volume of crater between minima 
%                       and highest boundary elevation.
%           Orientation: Structure of crater orientation parameters, 
%                   similar to OrientParams structure:
%               Major_Axis_Azimuth: Azimuthal direction of best-fitting 
%                   crater boundary ellipse.
%           Shape: Structure of crater shape parameters, similar to 
%                   ShapeParams structure:
%               Ellipticity_Index: Structure of crater ellipticity index,
%                       broken down by methodology.
%                   MaxDiameter: Ellipticity index of the crater using the 
%                       maximum diameter (from Grosse et al., 2012).
%                   BFEllipse: Ellipticity index  of the crater using the 
%                       best-fitting ellipse major axis.
%               Irregularity_Index: Structure of crater irregularity index,
%                       broken down by methodology.
%                   MaxDiameter: Irregularity index of the crater using the 
%                       maximum diameter (from Grosse et al., 2012).
%                   BFEllipse: Irregularity index  of the crater using the 
%                       best-fitting ellipse major axis.
%               Geometry: Structure of crater geometry values:
%                   CraterDepth_CraterWidth: Ratio between crater depth and 
%                       crater width.
%                   CraterWidth_BasalWidth: Ratio between crater width and 
%                       basal width.
%                   CraterDepth_BasalHeight: Ratio between crater depth and 
%                       edifice height.
%           Slope: Structure of crater shape parameters, similar to 
%                   SlopeParams structure:
%               Mean: Mean slope within crater.
%               Median: Median slope within crater.
%               Contour: Cell array of contour statistics within the crater:
%                   Elevations: Contour elevations.
%                   Mean: Mean slope along each contour.
%                   Median: Median slope along each contour.
%                   Min: Min slope along each contour.
%                   Max: Max slope along each contour.
%                   Std: Slope standard deviation along each contour.
%                   Maximum_Mean: Maximum mean slope of all contours.
%                   Maximum_Median: Maximum median slope of all contours.
%               Crater_Surfaces: Structure of filled-in crater surfaces 
%                   derived using the Natural, IDW, and Kriging methods 
%                   (based on user input). Given as GRIDobjs.
%
%       PeakParams: 
%           Full_Edifice: Structure of identified peaks for the entire 
%                   edifice, broken down by peak identification type.
%               Contour: Structure of peaks identified by the contour
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%               Local: Structure of peaks identified by the locality
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%           Lower_Flank: Structure of identified peaks for the lower 
%                   edifice flank (region between the lowest and highest 
%                   edifice boundary), broken down by peak identification 
%                   type.
%               Contour: Structure of peaks identified by the contour
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%               Local: Structure of peaks identified by the locality
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%           Main_Flank: Structure of identified peaks for the main 
%                   edifice flank (region between the highest edifice  
%                   boundar and summit), broken down by peak identification 
%                   type.
%               Contour: Structure of peaks identified by the contour
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%               Local: Structure of peaks identified by the locality
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%           Summit: Structure of identified peaks for the edifice summit, 
%                   broken down by peak identification type.
%               Contour: Structure of peaks identified by the contour
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.
%               Local: Structure of peaks identified by the locality
%                       method.
%                   Count: Number of peaks.
%                   XYZ: Peak X-Y-Z coordinates.

%% Unwrap Package
pack = MorVolc_DefaultVals(pack);
if pack.verbose > 0
    disp(sprintf('USING MORVOLC VERSION %s',pack.version))
    disp('Unpacking Input and Importing Data...')
end

MV_Res = [];
MV_Res.GeneralParams.inputs = pack;

dx = pack.dx;
summitRegion = pack.summitRegion;

MV_Res.GeneralParams.Version = pack.version;
MV_Res.GeneralParams.StartTime = datetime('now');

MorVolc_CheckIters(pack.contIter,pack.craterXY,pack.craterContIter,pack.peakContIter)
minEllipsePoints = 5;
ellIndMax = 3;
irrIndMax = 2;
irrIndMin = 1;

%% Import Shapefiles
[boundaryXYZ,craterXYZ,maskXY] = Import_Shapefiles(pack.boundaryXY,pack.craterXY,pack.maskXY,pack.verbose);

%% Load Tif file, isolate volcano
[DEM0,x0,y0,Z0,DEM,x,y,Z] = Import_DEM(pack.tifFile,dx,0,boundaryXYZ,craterXYZ,maskXY,pack.verbose);
[X,Y] = meshgrid(x,y);
[X0,Y0] = meshgrid(x0,y0);

%% Update Z values in boundary
bf = scatteredInterpolant(X0(:),Y0(:),Z0(:),'linear','linear');
bz = bf(boundaryXYZ(:,1),boundaryXYZ(:,2));
boundaryXYZ = [boundaryXYZ,bz];

for i = 1:length(craterXYZ)
    if pack.correctCrater
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
        for j = 1:size(bb{1}(:,1),1)
            txyz = [txyz;x(bb{1}(j,2)),y(bb{1}(j,1))];
        end 
    else
        txyz = craterXYZ{i}(:,1:2);
    end

    tz = interp2(X0,Y0,Z0,txyz(:,1),txyz(:,2));
    craterXYZ{i} = double([txyz,tz]);
end

MV_Res.GeographicParams.boundaryXYZ = boundaryXYZ;
MV_Res.GeographicParams.craterXYZ = craterXYZ;
MV_Res.GeographicParams.maskXY = maskXY;

MV_Res.GeographicParams.DEM0 = DEM0;
MV_Res.GeographicParams.DEM = DEM;

%% Create crater-masked slope grid
S = gradient8(DEM,'degree');
[SG,~,~] = GRIDobj2mat(S);
S_nocrater = SG;
for i = 1:length(craterXYZ) 
    txyz = craterXYZ{i}(:,1:2);

    p = inpolygon(X,Y,txyz(:,1),txyz(:,2));
    S_nocrater(p==1) = NaN;
end

MV_Res.GeographicParams.Slope = S;

%% Determine Contour Interval & Lower Flank Boundary
if pack.verbose > 0
    disp('Determing contours and lower flank boundary...')
end

if length(pack.contIter) == 1
    if pack.contIter < 0 && pack.contIter > -1
        trueCont = round(abs(pack.contIter)*range(Z(:)),2);
    else
        trueCont = round(pack.contIter,2);
    end

    % Determine contours
    tmpConts = round(nanmin(boundaryXYZ(:,3)),2):trueCont:round(nanmin(boundaryXYZ(:,3)),2)+trueCont;
    tmp = mod(tmpConts,trueCont);
    tmpI = find(tmp==0,1);
    if isempty(tmpI)
        tmpI = 1;
    end
    conts = [tmpConts(tmpI)-trueCont:trueCont:max(Z(:))+trueCont]';
else
    conts = zeros(length(pack.contIter),1);
    for i = 1:length(pack.contIter)
        if pack.contIter(i) < 0 && pack.contIter(i) > -1
            conts(i) = round(abs(pack.contIter(i))*range(Z(:)),2);
        else
            conts(i) = pack.contIter(i);
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
    LowFlankZ = nanmin(boundaryXYZ(:,3));
end

%% Generate basic topography metrics
if pack.verbose > 0
    disp('Calculating Topography Metrics...')
    disp('   Hypsometry...')
end

% Hypsome3try
[ZHyps_Areas,ZHyps_Vals] = HypsometryValue(DEM.Z,pack.hypsIter,1);
MV_Res.TopoParams.Hypsometry.Elevation_Hyps_Areas_Values = [ZHyps_Areas',ZHyps_Vals'];

% Slope
if pack.verbose > 0
    disp('   Slope...')
end
[Slope_Hyps_Areas,Slope_Hyps_Vals] = HypsometryValue(SG,pack.hypsIter,0);
MV_Res.TopoParams.Hypsometry.Slope_Hyps_Areas_Values = [Slope_Hyps_Areas',Slope_Hyps_Vals'];


% Slope Variance - Windows
if pack.verbose > 0
    disp('   Windowed Slope Variance...')
end
slopeVarianceVals = [];
for i = 1:length(pack.slopeVarianceWindows)
    if pack.verbose > 1
        disp(sprintf('      Window %d / %d',i,length(pack.slopeVarianceWindows)))
    end
    slVarI = round(pack.slopeVarianceWindows(i)/dx,0);
    if mod(slVarI,2) == 0
        slVarI = slVarI+1;
    end

    slVarGrid = CalculateSlopeVarianceWindow(S,slVarI);

    [slVarHypsAreas,slVarHypsVals] = HypsometryValue(slVarGrid.Z,pack.hypsIter,0);
    [~,slVarHypsNorm] = HypsometryValue(slVarGrid.Z,pack.hypsIter,1);
    
    slopeVarianceVals(i).SlopeVariance_Grid = slVarGrid;
    slopeVarianceVals(i).WindowSize = [slVarI,slVarI];
    slopeVarianceVals(i).ExpectedWindowRes = pack.slopeVarianceWindows(i);
    slopeVarianceVals(i).TrueWindowRes = slVarI*dx;
    slopeVarianceVals(i).Hypsometry_Areas = slVarHypsAreas;
    slopeVarianceVals(i).Hypsometry_Values = slVarHypsVals;
    slopeVarianceVals(i).Hypsometry_NormValues = slVarHypsNorm;
end

% Slope Variance - Contours
if pack.verbose > 0
    disp('   Contour Slope Variance...')
end
[slVarTable,slVarTitles] = CalculateSlopeVarianceElevation(DEM,S,conts);
MV_Res.TopoParams.SlopeVariance.Elevation.Values = slVarTable;
MV_Res.TopoParams.SlopeVariance.Elevation.Titles = slVarTitles;

MV_Res.TopoParams.SlopeVariance.Windows = slopeVarianceVals;
MV_Res.TopoParams.SlopeVariance.Total = nanstd(SG(:))./nanmean(SG(:));

% Curvature
if pack.verbose > 0
    disp('   Curvature...')
end
curve_prof = curvature(DEM,'profc');
curve_plan = curvature(DEM,'planc');

[CProf_Hyps_Areas,CProf_Hyps_Vals] = HypsometryValue(curve_prof.Z,pack.hypsIter*.1,0);
[CPlan_Hyps_Areas,CPlan_Hyps_Vals] = HypsometryValue(curve_plan.Z,pack.hypsIter*.1,0);

MV_Res.TopoParams.Curvature_Grids.Profile = curve_prof;
MV_Res.TopoParams.Curvature_Grids.Planform = curve_plan;

MV_Res.TopoParams.Hypsometry.ProfileC_Hyps_Areas_Values = [CProf_Hyps_Areas',CProf_Hyps_Vals'];
MV_Res.TopoParams.Hypsometry.PlanformC_Hyps_Areas_Values = [CPlan_Hyps_Areas',CPlan_Hyps_Vals'];

% Roughness
if pack.verbose > 0
    disp('   Roughness...')
end
roughnessVals = [];
for i = 1:length(pack.roughnessWindows)
    if pack.verbose > 1
        disp(sprintf('      Window %d / %d',i,length(pack.roughnessWindows)))
    end
    roughI = round(pack.roughnessWindows(i)/dx,0);
    if mod(roughI,2) == 0
        roughI = roughI+1;
    end
    roughGrid = roughness(DEM,pack.roughnessType,[roughI,roughI]);
    [roughHypsAreas,roughHypsVals] = HypsometryValue(roughGrid.Z,pack.hypsIter,0);
    [~,roughHypsNorm] = HypsometryValue(roughGrid.Z,pack.hypsIter,1);
    
    roughnessVals(i).Roughness_Grid = roughGrid;
    roughnessVals(i).WindowSize = [roughI,roughI];
    roughnessVals(i).ExpectedWindowRes = pack.roughnessWindows(i);
    roughnessVals(i).TrueWindowRes = roughI*dx;
    roughnessVals(i).Hypsometry_Areas = roughHypsAreas;
    roughnessVals(i).Hypsometry_Values = roughHypsVals;
    roughnessVals(i).Hypsometry_NormValues = roughHypsNorm;
end

MV_Res.TopoParams.Roughness = roughnessVals;

% Topographic, Geometric, and Volumetric Centers
%   Topographic Center (peak)
if pack.verbose > 0
    disp('   Landform Centers...')
end
[ii,jj] = find(Z==max(Z(:)),1);
TC_X = X(ii,jj);
TC_Y = Y(ii,jj);

%   Geometric Center (planform center of volcano)
%       The moments script calculates the volume-based mathmatical moments
%       of the elevation grid. Using 1 for all elevations within the
%       volcano removes the topographic weighting of the moments. This is
%       equivilant to taking the mean X and Y of all elevation coordinates.
tt = ~isnan(Z)*1;
[~,~,~,~,~,com,~,~] = Calculate_Moments(x,y,tt);
GC_X = com(1);
GC_Y = com(2);

%   Volumetric Center (volumetric center of mass, see Lerner et al. 2020 for
%       formal equation and description of use.
Z(isnan(Z)) = 0;
[~,~,~,~,~,com,~,~] = Calculate_Moments(x,y',Z);
Z(Z==0) = NaN;
VC_X = com(1);
VC_Y = com(2);

MV_Res.TopoParams.Centers.Topographic_XY = [TC_X,TC_Y];
MV_Res.TopoParams.Centers.Geometric_XY = [GC_X,GC_Y];
MV_Res.TopoParams.Centers.Volumetric_XY = [VC_X,VC_Y];

%% Determine Summit Area
if pack.verbose > 0
    disp('Determing summit area and contour parameters...')
end

if isempty(craterXYZ)
    useMaxZ = max(Z(:)); 
    if summitRegion==4 || summitRegion==5
        summitRegion = 1;
    end
else
    useMaxZ = max(Z(:));
    if summitRegion == 4
        for i = 1:length(craterXYZ)
            useMaxZ = min([useMaxZ,min(craterXYZ{i}(:,3))]);
        end
    elseif summitRegion == 5
        for i = 1:length(craterXYZ)
            tmpZ = Z;
            pp = inpolygon(X,Y,craterXYZ{i}(:,1),craterXYZ{i}(:,2));
            tmpZ(~pp) = NaN;
            useMaxZ = min([useMaxZ,nanmin(tmpZ(:))]);
        end
    end
end

[summitXYZ,summitCont] = MorVolc_GetSummitRegion(summitRegion,Z,SG,X,Y,LowFlankZ,conts,useMaxZ,pack.peakDiff,pack.verbose);

%% Mask grids and boundaries
if pack.verbose > 0
    disp('Applying Mask...')
end

if ~isempty(maskXY)
    for j = 1:length(maskXY)
        p = inpolygon(X,Y,maskXY{j}(:,1),maskXY{j}(:,2));
        tmpZ = Z;
        Z(p) = NaN;
        SG(p) = NaN;
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

MV_Res.GeographicParams.maskXY = maskXY;
MV_Res.GeographicParams.Lower_Flank_Contour = LowFlankZ;
MV_Res.GeographicParams.Lower_Flank_XY = LowFlankXYZ;

MV_Res.GeographicParams.Summit_Contour = summitCont;
MV_Res.GeographicParams.SummitXYZ = summitXYZ;
MV_Res.GeographicParams.boundaryXYZ = boundaryXYZ;
MV_Res.GeographicParams.craterXYZ = craterXYZ;

%% Calculate Contour Parameters
if pack.verbose > 0
    disp('Calculating Contour Parameters')
end

cont_Ellipses = cell(length(conts),1);
ellCenters = zeros(length(conts),2)*NaN;
contsForEI = cont_Ellipses;

cont_min_max_mean_median_Slope = zeros(length(conts),5)*NaN;

cont_min_roughnesses = zeros(length(conts),length(MV_Res.TopoParams.Roughness))*NaN;
cont_max_roughnesses = zeros(length(conts),length(MV_Res.TopoParams.Roughness))*NaN;
cont_mean_roughnesses = zeros(length(conts),length(MV_Res.TopoParams.Roughness))*NaN;
cont_median_roughnesses = zeros(length(conts),length(MV_Res.TopoParams.Roughness))*NaN;
cont_std_roughnesses = zeros(length(conts),length(MV_Res.TopoParams.Roughness))*NaN;

cont_min_sv = zeros(length(conts),length(MV_Res.TopoParams.SlopeVariance.Windows))*NaN;
cont_max_sv = zeros(length(conts),length(MV_Res.TopoParams.SlopeVariance.Windows))*NaN;
cont_mean_sv = zeros(length(conts),length(MV_Res.TopoParams.SlopeVariance.Windows))*NaN;
cont_median_sv = zeros(length(conts),length(MV_Res.TopoParams.SlopeVariance.Windows))*NaN;
cont_std_sv = zeros(length(conts),length(MV_Res.TopoParams.SlopeVariance.Windows))*NaN;

R_Grids = cell(length(MV_Res.TopoParams.Roughness),1);
for i = 1:length(R_Grids)
    [R,~,~] = GRIDobj2mat(MV_Res.TopoParams.Roughness(i).Roughness_Grid);
    R(isnan(S_nocrater)) = NaN;
    R_Grids{i} = R;
end

SV_Grids = cell(length(MV_Res.TopoParams.SlopeVariance.Windows),1);
for i = 1:length(SV_Grids)
    [SV,~,~] = GRIDobj2mat(MV_Res.TopoParams.SlopeVariance.Windows(i).SlopeVariance_Grid);
    SV(isnan(S_nocrater)) = NaN;
    SV_Grids{i} = SV;
end

if pack.verbose > 0
    disp('   Determining ellipses, slopes, roughnesses, and slope variance...')
end
for i = 1:length(conts)
    if pack.verbose > 1
        disp(sprintf('      %d / %d',i,length(conts)));
    end

    % Seperate contour
    ccXY = contourc(X(1,:),Y(:,1),Z,[1,1]*conts(i));
    ccXY_All = Convert_Contours(ccXY,0);
    ccXY_Sing = Convert_Contours(ccXY,1);

    if isempty(ccXY)
        continue;
    end

    % Get and fill in slopes of the contour
    ss = S_nocrater;
    ss(Z>conts(i)) = NaN;
    ss(Z<=conts(i-1)) = NaN;

    cont_min_max_mean_median_Slope(i,1) = nanmin(ss(:));
    cont_min_max_mean_median_Slope(i,2) = nanmax(ss(:));
    cont_min_max_mean_median_Slope(i,3) = nanmean(ss(:));
    cont_min_max_mean_median_Slope(i,4) = nanmedian(ss(:));
    cont_min_max_mean_median_Slope(i,5) = nanstd(ss(:));

    % Determine ellipse of closed contours
    if conts(i) >= LowFlankZ
        contsForEI{i} = ccXY_Sing;
        ce = FitEllipse_Unwrap(ccXY_Sing,[]);
        if ~isreal(ce.shortAxis)
            ce.shortAxis = NaN;
            ce.longAxis = NaN;
            ce.phi = NaN;
        else
            ellCenters(i,:) = [ce.x0,ce.y0];
        end
        cont_Ellipses{i} = ce;
    else
        contsForEI{i} = ccXY_All;
    end

    if i == 1
        continue;
    end

    % Loop through roughness, get and fill in contour values
    for j = 1:length(R_Grids)
        R = R_Grids{j};
        R(Z>conts(i)) = NaN;
        R(Z<=conts(i-1)) = NaN;
        
        cont_min_roughnesses(i,j) = nanmin(R(:));
        cont_max_roughnesses(i,j) = nanmax(R(:));
        cont_mean_roughnesses(i,j) = nanmean(R(:));
        cont_median_roughnesses(i,j) = nanmedian(R(:));
        cont_std_roughnesses(i,j) = nanstd(R(:));
    end

    % Loop through slope variance, get and fill in contour values
    for j = 1:length(SV_Grids)
        SV = SV_Grids{j};
        SV(Z>conts(i)) = NaN;
        SV(Z<=conts(i-1)) = NaN;
        
        cont_min_sv(i,j) = nanmin(SV(:));
        cont_max_sv(i,j) = nanmax(SV(:));
        cont_mean_sv(i,j) = nanmean(SV(:));
        cont_median_sv(i,j) = nanmedian(SV(:));
        cont_std_sv(i,j) = nanstd(SV(:));
    end
end

% Try to fill in ellipses of non-closed contours
if pack.verbose > 0
    disp('   Correcting ellipses...')
end
meanXY = nanmean(ellCenters,1);
for i = 1:length(conts)
    if ~isempty(cont_Ellipses{i})
        continue;
    end

    contPoints = [];
    consecutivePointsHit = 0;
    for j = 1:length(contsForEI{i})
        if size(contsForEI{i}{j},1) >= minEllipsePoints
            consecutivePointsHit = 1;
        end
        contPoints = [contPoints;contsForEI{i}{j}];
    end

    if ~consecutivePointsHit
        continue;
    end

    ce = FitEllipse_Unwrap(contPoints,meanXY);
    if ~isreal(ce.shortAxis)
        ce.shortAxis = NaN;
        ce.longAxis = NaN;
        ce.phi = NaN;
    end
    cont_Ellipses{i} = ce;
end

% Calculate EI and II
if pack.verbose > 0
    disp('   Calculating ellipticity and irregularity indices...')
end

ei_ii_md = zeros(length(conts),2)*NaN;
ei_ii_bfe = ei_ii_md;
ee = zeros(length(conts),1)*NaN;
azS = ee;
for i = 1:length(conts)
    if ~isempty(cont_Ellipses{i})
        ee(i) = cont_Ellipses{i}.shortAxis/cont_Ellipses{i}.longAxis;
        azS(i) = cont_Ellipses{i}.phi;
    end

    if isempty(contsForEI{i}) || (isempty(cont_Ellipses{i}) && conts(i)<LowFlankZ) || (pack.ignoreSummitIndex==0 && conts(i) > min(summitXYZ(:,3))) 
        continue;
    end

    [ei_md,ii_md,~,~] = MorVolc_EI_II(contsForEI{i},cont_Ellipses{i},conts(i)>=LowFlankZ,0);
    [ei_bfe,ii_bfe,~,~] = MorVolc_EI_II(contsForEI{i},cont_Ellipses{i},conts(i)>=LowFlankZ,1);

    ei_ii_md(i,1) = ei_md;
    ei_ii_md(i,2) = ii_md;
    ei_ii_bfe(i,1) = ei_bfe;
    ei_ii_bfe(i,2) = ii_bfe;
end

%% Collect Size Parameters
if pack.verbose > 0
    disp('Collecting Size Parameters...')
end
% Basal Area & Width
[basalArea,basalWidth] = MorVolc_Calculate_Area_Width(boundaryXYZ(:,1),boundaryXYZ(:,2));

MV_Res.SizeParams.Basal.Area = basalArea;
MV_Res.SizeParams.Basal.Width = basalWidth;

% Basal Axes
be = FitEllipse_Unwrap(boundaryXYZ(:,1:2),[]);
MV_Res.SizeParams.Basal.Major_Axis = be.longAxis;
MV_Res.SizeParams.Basal.Minor_Axis = be.shortAxis;
MV_Res.SizeParams.Basal.Axis_Ellipticity = be.shortAxis/be.longAxis;
MV_Res.SizeParams.Basal.Ellipse = be;

% Summit Area & Width
[summitArea,summitWidth] = MorVolc_Calculate_Area_Width(summitXYZ(:,1),summitXYZ(:,2));

MV_Res.SizeParams.Summit.Area = summitArea;
MV_Res.SizeParams.Summit.Width = summitWidth;

% Summit Axes
se = FitEllipse_Unwrap(summitXYZ(:,1:2));
MV_Res.SizeParams.Summit.Major_Axis = se.longAxis;
MV_Res.SizeParams.Summit.Minor_Axis = se.shortAxis;
MV_Res.SizeParams.Summit.Axis_Ellipticity = se.shortAxis/se.longAxis;
MV_Res.SizeParams.Summit.Ellipse = se;

% Heights and Volumes
[~,~,surfDEMs,maxHeight,maxVol,volcHeights,volcVols,anyHeights] = MorVolc_CalculateHeightVols(X,Y,Z,boundaryXYZ,min(boundaryXYZ(:,3)),pack.interpSurfaces);

MV_Res.SizeParams.Heights.Surface_to_Peak = volcHeights;
MV_Res.SizeParams.Heights.Maximum_From_Boundary = maxHeight;
MV_Res.SizeParams.Heights.Maximum_From_BasalSurface = anyHeights;

MV_Res.GeographicParams.Basal_Surfaces = surfDEMs;
% MV_Res.GeographicParams.Basal_Surface_XY = surfXYs;
% MV_Res.GeographicParams.Basal_Surface_Z = surfZs;

MV_Res.SizeParams.Volumes.Total = volcVols;
MV_Res.SizeParams.Volumes.Maximum = maxVol;

% Eroded Volume
[cont_area_convHullArea,volumeDiff,ConvZ] = CalculateErodedVolume(X,Y,Z,boundaryXYZ,trueCont,min(boundaryXYZ(:,3)),1,pack.verbose);
MV_Res.SizeParams.Volumes.Minimum_Eroded = volumeDiff;
MV_Res.SizeParams.ReconstructedTopo.Convex_Hull_Areas = cont_area_convHullArea;
MV_Res.SizeParams.ReconstructedTopo.Convex_Hull_Interpolated_Surface = ConvZ;

%% Collect Orientation Parameters
if pack.verbose > 0
    disp('Collecting Orientation Parameters...')
end
MV_Res.OrientParams.Basal_Major_Axis_Azimuth = be.phi;
MV_Res.OrientParams.Summit_Major_Axis_Azimuth = se.phi;
MV_Res.OrientParams.Contour_Elevation_Major_Axis_Azimuths = [conts,azS];

%% Collect Shape Parameters
if pack.verbose > 0
    disp('Collecting Shape Parameters...')
end
MV_Res.ShapeParams.Contour.Elevations = conts;
MV_Res.ShapeParams.Contour.Ellipses = cont_Ellipses;
MV_Res.ShapeParams.Contour.Axis_Ellipticity = ee;

ei_ii_md(ei_ii_md(:,1)>ellIndMax,1) = ellIndMax;
ei_ii_md(ei_ii_md(:,2)>irrIndMax,2) = irrIndMax;
ei_ii_md(ei_ii_md(:,2)<irrIndMin,2) = NaN;

MV_Res.ShapeParams.Contour.Ellipticity_Indices.MaxDiameter = ei_ii_md(:,1);
MV_Res.ShapeParams.Contour.Irregularity_Indices.MaxDiameter = ei_ii_md(:,2);

MV_Res.ShapeParams.Mean_Shapes.Ellipticity_Index.MaxDiameter = nanmean(ei_ii_md(:,1));
MV_Res.ShapeParams.Mean_Shapes.Irregularity_Index.MaxDiameter = nanmean(ei_ii_md(:,2));

ei_ii_bfe(ei_ii_bfe(:,1)>ellIndMax,1) = ellIndMax;
ei_ii_bfe(ei_ii_bfe(:,2)>irrIndMax,2) = irrIndMax;
ei_ii_bfe(ei_ii_bfe(:,2)<irrIndMin,2) = NaN;

MV_Res.ShapeParams.Contour.Ellipticity_Indices.BFEllipse = ei_ii_bfe(:,1);
MV_Res.ShapeParams.Contour.Irregularity_Indices.BFEllipse = ei_ii_bfe(:,2);

MV_Res.ShapeParams.Mean_Shapes.Ellipticity_Index.BFEllipse = nanmean(ei_ii_bfe(:,1));
MV_Res.ShapeParams.Mean_Shapes.Irregularity_Index.BFEllipse = nanmean(ei_ii_bfe(:,2));

MV_Res.ShapeParams.Mean_Shapes.Ellipse_Ellipticity = nanmean(ee);

MV_Res.ShapeParams.Geometry.Height_BasalWidth.Natural = MV_Res.SizeParams.Heights.Surface_to_Peak.Natural/MV_Res.SizeParams.Basal.Width;
MV_Res.ShapeParams.Geometry.Height_BasalWidth.IDW = MV_Res.SizeParams.Heights.Surface_to_Peak.IDW/MV_Res.SizeParams.Basal.Width;
MV_Res.ShapeParams.Geometry.Height_BasalWidth.Kriging = MV_Res.SizeParams.Heights.Surface_to_Peak.Kriging/MV_Res.SizeParams.Basal.Width;
MV_Res.ShapeParams.Geometry.SummitWidth_BasalWidth = MV_Res.SizeParams.Summit.Width/MV_Res.SizeParams.Basal.Width;

Znan = isnan(Z);
Z(isnan(Z)) = 0;
[~,~,~,~,~,~,sk,ku] = Calculate_Moments(x,y',Z);
Z(Znan) = NaN;
MV_Res.ShapeParams.Geometry.Skewness = sk;
MV_Res.ShapeParams.Geometry.Kurtosis = ku;

%% Collect Slope Parameters
if pack.verbose > 0
    disp('Collecting Slope Parameters...')
end
tmpS = S_nocrater;
tmpS(isnan(Z)) = NaN;
MV_Res.SlopeParams.Full_Edifice.Mean = nanmean(tmpS(:));
MV_Res.SlopeParams.Full_Edifice.Median = nanmedian(tmpS(:));
MV_Res.SlopeParams.Full_Edifice.Std = nanstd(tmpS(:));

pF = inpolygon(X,Y,LowFlankXYZ(:,1),LowFlankXYZ(:,2));
pS = inpolygon(X,Y,summitXYZ(:,1),summitXYZ(:,2));

tmpS = S_nocrater;
tmpS(pF) = NaN;
MV_Res.SlopeParams.Lower_Flank.Mean = nanmean(tmpS(:));
MV_Res.SlopeParams.Lower_Flank.Median = nanmedian(tmpS(:));
MV_Res.SlopeParams.Lower_Flank.Std = nanstd(tmpS(:));

tmpS = S_nocrater;
tmpS(~pF) = NaN;
tmpS(pS) = NaN;
MV_Res.SlopeParams.Main_Flank.Mean = nanmean(tmpS(:));
MV_Res.SlopeParams.Main_Flank.Median = nanmedian(tmpS(:));
MV_Res.SlopeParams.Main_Flank.Std = nanstd(tmpS(:));

tmpS = S_nocrater;
tmpS(~pS) = NaN;
MV_Res.SlopeParams.Summit.Mean = nanmean(tmpS(:));
MV_Res.SlopeParams.Summit.Median = nanmedian(tmpS(:));
MV_Res.SlopeParams.Summit.Std = nanstd(tmpS(:));

MV_Res.SlopeParams.Contour.Elevations = conts;
MV_Res.SlopeParams.Contour.Mean = cont_min_max_mean_median_Slope(:,3);
MV_Res.SlopeParams.Contour.Median = cont_min_max_mean_median_Slope(:,4);
MV_Res.SlopeParams.Contour.Min = cont_min_max_mean_median_Slope(:,1);
MV_Res.SlopeParams.Contour.Max = cont_min_max_mean_median_Slope(:,2);
MV_Res.SlopeParams.Contour.Std = cont_min_max_mean_median_Slope(:,5);
MV_Res.SlopeParams.Contour.Maximum_Mean = max(cont_min_max_mean_median_Slope(:,3));
MV_Res.SlopeParams.Contour.Maximum_Median = max(cont_min_max_mean_median_Slope(:,4));

hf_i = find(cont_min_max_mean_median_Slope(:,2)==max(cont_min_max_mean_median_Slope(:,2)));
el_cont = conts(hf_i);
hf_cont= el_cont-min(boundaryXYZ(:,3));
totHeight = max(Z(:)) - min(boundaryXYZ(:,3));
heightFrac = hf_cont/totHeight;

MV_Res.SlopeParams.Maximum_Slope.Elevation = el_cont;
MV_Res.SlopeParams.Maximum_Slope.Height = hf_cont;
MV_Res.SlopeParams.Maximum_Slope.Height_Fraction = heightFrac;
MV_Res.SlopeParams.Maximum_Slope.Slope = max(cont_min_max_mean_median_Slope(:,2));

%% Collect Roughness Parameters
if pack.verbose > 0
    disp('Collecting Roughness Parameters...')
end
tmp = zeros(1,length(R_Grids))*NaN;
MV_Res.RoughnessParams.Full_Edifice.Mean = tmp;
MV_Res.RoughnessParams.Full_Edifice.Median = tmp;
MV_Res.RoughnessParams.Full_Edifice.Std = tmp;

MV_Res.RoughnessParams.Lower_Flank.Mean = tmp;
MV_Res.RoughnessParams.Lower_Flank.Median = tmp;
MV_Res.RoughnessParams.Lower_Flank.Std = tmp;

MV_Res.RoughnessParams.Main_Flank.Mean = tmp;
MV_Res.RoughnessParams.Main_Flank.Median = tmp;
MV_Res.RoughnessParams.Main_Flank.Std = tmp;

MV_Res.RoughnessParams.Summit.Mean = tmp;
MV_Res.RoughnessParams.Summit.Median = tmp;
MV_Res.RoughnessParams.Summit.Std = tmp;

for i = 1:length(R_Grids)
    tmpR = R_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(isnan(Z)) = NaN;
    MV_Res.RoughnessParams.Full_Edifice.Mean(i) = nanmean(tmpR(:));
    MV_Res.RoughnessParams.Full_Edifice.Median(i) = nanmedian(tmpR(:));
    MV_Res.RoughnessParams.Full_Edifice.Std(i) = nanstd(tmpR(:));

    tmpR = R_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(pF) = NaN;
    MV_Res.RoughnessParams.Lower_Flank.Mean(i) = nanmean(tmpR(:));
    MV_Res.RoughnessParams.Lower_Flank.Median(i) = nanmedian(tmpR(:));
    MV_Res.RoughnessParams.Lower_Flank.Std(i) = nanstd(tmpR(:));

    tmpR = R_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(~pS) = NaN;
    MV_Res.RoughnessParams.Main_Flank.Mean(i) = nanmean(tmpR(:));
    MV_Res.RoughnessParams.Main_Flank.Median(i) = nanmedian(tmpR(:));
    MV_Res.RoughnessParams.Main_Flank.Std(i) = nanstd(tmpR(:));

    tmpR = R_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(~pF) = NaN;
    tmpR(pS) = NaN;
    MV_Res.RoughnessParams.Summit.Mean(i) = nanmean(tmpR(:));
    MV_Res.RoughnessParams.Summit.Median(i) = nanmedian(tmpR(:));
    MV_Res.RoughnessParams.Summit.Std(i) = nanstd(tmpR(:));
end

MV_Res.RoughnessParams.Contour.Elevations = conts;
MV_Res.RoughnessParams.Contour.Mean = cont_mean_roughnesses;
MV_Res.RoughnessParams.Contour.Median = cont_median_roughnesses;
MV_Res.RoughnessParams.Contour.Min = cont_min_roughnesses;
MV_Res.RoughnessParams.Contour.Max = cont_max_roughnesses;
MV_Res.RoughnessParams.Contour.Std = cont_std_roughnesses;
MV_Res.RoughnessParams.Contour.Maximum_Mean = nanmax(cont_mean_roughnesses,[],1);
MV_Res.RoughnessParams.Contour.Maximum_Median = nanmax(cont_median_roughnesses,[],1);

%% Collect Slope Variance Parameters
if pack.verbose > 0
    disp('Collecting Windowed Slope Variance Parameters...')
end
% tmp = zeros(1,length(R_Grids))*NaN;
tmp = zeros(1,length(SV_Grids))*NaN;
MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Mean = tmp;
MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Median = tmp;
MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Std = tmp;

MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Mean = tmp;
MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Median = tmp;
MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Std = tmp;

MV_Res.WindowedSlopeVarianceParams.Main_Flank.Mean = tmp;
MV_Res.WindowedSlopeVarianceParams.Main_Flank.Median = tmp;
MV_Res.WindowedSlopeVarianceParams.Main_Flank.Std = tmp;

MV_Res.WindowedSlopeVarianceParams.Summit.Mean = tmp;
MV_Res.WindowedSlopeVarianceParams.Summit.Median = tmp;
MV_Res.WindowedSlopeVarianceParams.Summit.Std = tmp;

for i = 1:length(SV_Grids)
    tmpR = SV_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(isnan(Z)) = NaN;
    MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Mean(i) = nanmean(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Median(i) = nanmedian(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Full_Edifice.Std(i) = nanstd(tmpR(:));

    % tmpR = R_Grids{i};
    tmpR = SV_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(pF) = NaN;
    MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Mean(i) = nanmean(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Median(i) = nanmedian(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Lower_Flank.Std(i) = nanstd(tmpR(:));

    % tmpR = R_Grids{i};
    tmpR = SV_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(~pS) = NaN;
    MV_Res.WindowedSlopeVarianceParams.Main_Flank.Mean(i) = nanmean(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Main_Flank.Median(i) = nanmedian(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Main_Flank.Std(i) = nanstd(tmpR(:));

    % tmpR = R_Grids{i};
    tmpR = SV_Grids{i};
    tmpR(isnan(S_nocrater)) = NaN;
    tmpR(~pF) = NaN;
    tmpR(pS) = NaN;
    MV_Res.WindowedSlopeVarianceParams.Summit.Mean(i) = nanmean(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Summit.Median(i) = nanmedian(tmpR(:));
    MV_Res.WindowedSlopeVarianceParams.Summit.Std(i) = nanstd(tmpR(:));
end

MV_Res.WindowedSlopeVarianceParams.Contour.Elevations = conts;
MV_Res.WindowedSlopeVarianceParams.Contour.Mean = cont_mean_sv;
MV_Res.WindowedSlopeVarianceParams.Contour.Median = cont_median_sv;
MV_Res.WindowedSlopeVarianceParams.Contour.Min = cont_min_sv;
MV_Res.WindowedSlopeVarianceParams.Contour.Max = cont_max_sv;
MV_Res.WindowedSlopeVarianceParams.Contour.Std = cont_std_sv;
MV_Res.WindowedSlopeVarianceParams.Contour.Maximum_Mean = nanmax(cont_mean_sv,[],1);
MV_Res.WindowedSlopeVarianceParams.Contour.Maximum_Median = nanmax(cont_median_sv,[],1);

%% Collect Crater Parameters
if pack.verbose > 0
    disp('Collecting Crater Parameters...')
end

if isempty(craterXYZ)
    MV_Res.CraterParams = [];
    % MV_Res.GeneralParams.Crater_Interp = [];
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
    % craterSurfXY = craterEllipse;
    % craterSurfZ = craterEllipse;
    craterDEMs = craterEllipse;
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
        [ca,cw] = MorVolc_Calculate_Area_Width(cxyz(:,1),cxyz(:,2));
        
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
        [cEi_MaxDiam,cIi_MaxDiam,~,~] = MorVolc_EI_II(cxyz,cEl,1,0);
        [cEi_BFEllipse,cIi_BFEllipse,~,~] = MorVolc_EI_II(cxyz,cEl,1,1);
        
        craterElipt_MaxDiam(i) = cEi_MaxDiam;
        craterIreg_MaxDiam(i) = cIi_MaxDiam;
        craterElipt_BFEllipse(i) = cEi_BFEllipse;
        craterIreg_BFEllipse(i) = cIi_BFEllipse;
        
        % Depth & Volume
        tmpZ = Z;
        p = inpolygon(X,Y,cxyz(:,1),cxyz(:,2));
        tmpZ(~p) = NaN;
        tmpCxyz = cxyz;
        [cSurfXYs,cSurfZs,cDEMs,mH,mV,vHs,vVs] = MorVolc_CalculateHeightVols_Crater(X,Y,tmpZ,tmpCxyz,max(tmpCxyz(:,3)),pack.interpSurfaces);

        craterMaxDepth(i) = mH;
        craterMaxVol(i) = mV;
        craterDepth{i} = vHs;
        craterVol{i} = vVs;
        % craterSurfXY{i} = cSurfXYs;
        % craterSurfZ{i} = cSurfZs;
        craterDEMs{i} = cDEMs;
  
        % Slopes
        tmpS = SG;
        tmpS(~p) = NaN;
        craterMeanSl(i) = nanmean(tmpS(:));
        craterMedianSl(i) = nanmedian(tmpS(:));

        tmpZ = Z;
        tmpZ(~p) = NaN;
        
        % Contour Stats
        if pack.craterContIter < 0 && pack.craterContIter > -1
            trueCraterCont = round(abs(pack.craterContIter)*range(tmpZ(:)));
        else
            trueCraterCont = round(pack.craterContIter);
        end
        
        tmpConts = round(min(tmpZ(:))):round(min(tmpZ(:)))+trueCraterCont;
        tmp = mod(tmpConts,trueCraterCont);
        tmpI = find(tmp==0,1);
        craterConts = [tmpConts(tmpI)-trueCraterCont:trueCraterCont:max(tmpZ(:))]';
        
        
        tmpStats.Elevations = craterConts;
        tmpCMeanSl = zeros(size(tmpStats.Elevations));
        tmpCMedianSl = zeros(size(tmpStats.Elevations));
        tmpCMinSl = zeros(size(tmpStats.Elevations));
        tmpCMaxSl = zeros(size(tmpStats.Elevations));
        tmpCStdSl = zeros(size(tmpStats.Elevations));
        
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
            tmpCStdSl(j) = nanstd(ss(:));
        end
        
        tmpStats.Mean = tmpCMeanSl;
        tmpStats.Median = tmpCMedianSl;
        tmpStats.Min = tmpCMinSl;
        tmpStats.Max = tmpCMaxSl;
        tmpStats.Std = tmpCStdSl;
        tmpStats.Maximum_Mean = nanmax(tmpCMeanSl);
        tmpStats.Maximum_Median = nanmax(tmpCMedianSl);
        
        craterContStats{i} = tmpStats;
        
        % Other Crater Stats
        craterDepth_Width{i}.Natural = craterDepth{i}.Natural/craterWidth(i);
        craterDepth_Width{i}.IDW = craterDepth{i}.IDW/craterWidth(i);
        craterDepth_Width{i}.Kriging = craterDepth{i}.Kriging/craterWidth(i);
        craterWidth_BasalWidth(i) = craterWidth(i)/basalWidth;
        craterDepth_BasalHeight{i}.Natural = craterDepth{i}.Natural/volcHeights.Natural;
        craterDepth_BasalHeight{i}.IDW = craterDepth{i}.IDW/volcHeights.IDW;
        craterDepth_BasalHeight{i}.Kriging = craterDepth{i}.Kriging/volcHeights.Kriging;
    end
    
    MV_Res.CraterParams.Size.Area = craterArea;
    MV_Res.CraterParams.Size.Width = craterWidth;
    MV_Res.CraterParams.Size.Major_Axis = craterMajor;
    MV_Res.CraterParams.Size.Minor_Axis = craterMinor;
    MV_Res.CraterParams.Size.Axis_Ellipticity = craterAxisEllipticity;
    MV_Res.CraterParams.Size.Ellipse = craterEllipse;
    MV_Res.CraterParams.Size.Depth.Total = craterDepth;
    MV_Res.CraterParams.Size.Depth.Maximum = craterMaxDepth;
    MV_Res.CraterParams.Size.Volume.Total = craterVol;
    MV_Res.CraterParams.Size.Volume.Maximum = craterMaxVol;

    MV_Res.CraterParams.Orientation.Major_Axis_Azimuth = craterAz;

    MV_Res.CraterParams.Shape.Ellipticity_Index.MaxDiameter = craterElipt_MaxDiam;
    MV_Res.CraterParams.Shape.Irregularity_Index.MaxDiameter = craterIreg_MaxDiam;
    MV_Res.CraterParams.Shape.Ellipticity_Index.BFEllipse = craterElipt_BFEllipse;
    MV_Res.CraterParams.Shape.Irregularity_Index.BFEllipse = craterIreg_BFEllipse;

    MV_Res.CraterParams.Shape.Geometry.CraterDepth_CraterWidth = craterDepth_Width;
    MV_Res.CraterParams.Shape.Geometry.CraterWidth_BasalWidth = craterWidth_BasalWidth;
    MV_Res.CraterParams.Shape.Geometry.CraterDepth_BasalHeight = craterDepth_BasalHeight;
    
    MV_Res.CraterParams.Slope.Mean = craterMeanSl;
    MV_Res.CraterParams.Slope.Median = craterMedianSl;
    MV_Res.CraterParams.Slope.Contour = craterContStats;
    
    MV_Res.CraterParams.Crater_Surfaces = craterDEMs;
end

%% Collect Peak Parameters
if pack.verbose > 0
    disp('Collecting Peak Parameters...')
end

if pack.peakContIter < 0 && pack.peakContIter > -1
    truePeakCont = abs(pack.peakContIter)*maxHeight;
else
    truePeakCont = pack.peakContIter;
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

allPeakXYZ_cont = zeros(size(allPeakIJ_cont,1),3);
for i = 1:size(allPeakIJ_cont,1)
    allPeakXYZ_cont(i,:) = [x(allPeakIJ_cont(i,2)),y(allPeakIJ_cont(i,1)),Z(allPeakIJ_cont(i,1),allPeakIJ_cont(i,2))];
end

[sinp,son] = inpolygon(allPeakXYZ_cont(:,1),allPeakXYZ_cont(:,2),summitXYZ(:,1),summitXYZ(:,2));
[lin,lon] = inpolygon(allPeakXYZ_cont(:,1),allPeakXYZ_cont(:,2),LowFlankXYZ(:,1),LowFlankXYZ(:,2));

MV_Res.PeakParams.Full_Edifice = [];
MV_Res.PeakParams.Lower_Flank.Contour.Count = size(allPeakXYZ_cont,1)-sum((lin+lon)>0);
MV_Res.PeakParams.Lower_Flank.Contour.XYZ = allPeakXYZ_cont((lin+lon)==0,:);
MV_Res.PeakParams.Main_Flank.Contour.Count = sum((lin+lon)>0)-sum((sinp+son)>0);
MV_Res.PeakParams.Main_Flank.Contour.XYZ = allPeakXYZ_cont(((sinp+son)==0).*((lin+lon)>0)>0,:);
MV_Res.PeakParams.Summit.Contour.Count = sum((sinp+son)>0);
MV_Res.PeakParams.Summit.Contour.XYZ = allPeakXYZ_cont((sinp+son)>0,:);

MV_Res.PeakParams.Full_Edifice.Contour.Count = MV_Res.PeakParams.Lower_Flank.Contour.Count +...
    MV_Res.PeakParams.Main_Flank.Contour.Count + MV_Res.PeakParams.Summit.Contour.Count;
MV_Res.PeakParams.Full_Edifice.Contour.XYZ = [MV_Res.PeakParams.Lower_Flank.Contour.XYZ;...
    MV_Res.PeakParams.Main_Flank.Contour.XYZ;MV_Res.PeakParams.Summit.Contour.XYZ];

% Local Maxima Peak Count
locMax1 = islocalmax(Z,1);
locMax2 = islocalmax(Z,2);
locMax3 = locMax1.*locMax2;

lM3_summit = locMax3;
[sinp,son] = inpolygon(X,Y,summitXYZ(:,1),summitXYZ(:,2));
lM3_summit((sinp+son)==0) = NaN;

lM3_mFlank = locMax3;
[lin,lon] = inpolygon(X,Y,LowFlankXYZ(:,1),LowFlankXYZ(:,2));
lM3_mFlank((lin+lon)==0) = NaN;
lM3_mFlank((sinp+son)>0) = NaN;

lM3_lFlank = locMax3;
lM3_lFlank((lin+lon)>0) = NaN;

MV_Res.PeakParams.Summit.Local.Count = nansum(lM3_summit(:));
[ii,jj] = find(lM3_summit>0);
MV_Res.PeakParams.Summit.Local.XYZ = [];
for i = 1:length(ii)
    MV_Res.PeakParams.Summit.Local.XYZ = [MV_Res.PeakParams.Summit.Local.XYZ;...
        X(ii(i),jj(i)),Y(ii(i),jj(i)),Z(ii(i),jj(i))];
end

MV_Res.PeakParams.Main_Flank.Local.Count = nansum(lM3_mFlank(:));
[ii,jj] = find(lM3_mFlank>0);
MV_Res.PeakParams.Main_Flank.Local.XYZ = [];
for i = 1:length(ii)
    MV_Res.PeakParams.Main_Flank.Local.XYZ = [MV_Res.PeakParams.Main_Flank.Local.XYZ;...
        X(ii(i),jj(i)),Y(ii(i),jj(i)),Z(ii(i),jj(i))];
end

MV_Res.PeakParams.Lower_Flank.Local.Count = nansum(lM3_lFlank(:));
[ii,jj] = find(lM3_lFlank>0);
MV_Res.PeakParams.Lower_Flank.Local.XYZ = [];
for i = 1:length(ii)
    MV_Res.PeakParams.Lower_Flank.Local.XYZ = [MV_Res.PeakParams.Lower_Flank.Local.XYZ;...
        X(ii(i),jj(i)),Y(ii(i),jj(i)),Z(ii(i),jj(i))];
end

MV_Res.PeakParams.Full_Edifice.Local.Count = MV_Res.PeakParams.Lower_Flank.Local.Count +...
    MV_Res.PeakParams.Main_Flank.Local.Count + MV_Res.PeakParams.Summit.Local.Count;
MV_Res.PeakParams.Full_Edifice.Local.XYZ = [MV_Res.PeakParams.Lower_Flank.Local.XYZ;...
    MV_Res.PeakParams.Main_Flank.Local.XYZ;MV_Res.PeakParams.Summit.Local.XYZ];

%% Save Results
MV_Res.GeneralParams.EndTime = datetime('now');
if ~isempty(pack.saveResFolder)
    if pack.verbose > 0
        disp('Saving Results...')
    end
    save([pack.saveResFolder,pack.figPrefix,'MorVolc_Results.mat'],'MV_Res')
end

if ~isempty(pack.xlsFile) && pack.xlsFile && ~isempty(pack.saveResFolder)
    MorVolc_SaveXLS([pack.saveResFolder,pack.figPrefix,'MorVolc_Results.xlsx'],MV_Res);
end

%% Save Inputs
if pack.saveInputs && ~isempty(pack.saveResFolder)
    if pack.verbose > 0
        disp('Writing Input Text File...')
    end

    tmpTif = pack.tifFile;
    tmpMask = pack.maskXY;
    tmpCrat = pack.craterXY;
    tmpXls = pack.xlsFile;

    if ~ischar(tmpTif)
        pack.tifFile = 'GRIDobj';
    end

    if isempty(tmpMask)
        pack.maskXY = '[]';
    elseif iscell(tmpMask)
        pack.maskXY = 'Cell array';
    elseif ~ischar(tmpMask)
        pack.maskXY = 'Array';
    end

    if isempty(tmpCrat)
        pack.craterXY = '[]';
    elseif iscell(tmpCrat)
        pack.craterXY = 'Cell array';
    elseif ~ischar(tmpCrat)
        pack.craterXY = 'Array';
    end

    if isempty(tmpXls)
        pack.xlsFile = '[]';
    end

    Export_Inputs(pack,'MorVolc');

    pack.tifFile = tmpTif;
    pack.maskXY = tmpMask;
    pack.craterXY = tmpCrat;
    pack.xlsFile = tmpXls;
end

%% Plot Results
if pack.plotResults
    if pack.verbose > 0
        disp('Plotting Results...')
    end
    MorVolc_Plots(MV_Res);
end

%% Zipping and Folder Deletion
if pack.zipFiles > 0
    if pack.verbose > 0   
        disp('Zipping Results...')
    end

    origPath = pwd;
    sameFol = (~isempty(pack.saveFigFolder) && ~isempty(pack.saveResFolder) && strcmp(pack.saveFigFolder,pack.saveResFolder));

    if sameFol || (~isempty(pack.saveFigFolder) && (pack.zipFiles == 1 || pack.zipFiles == 3))
        if pack.verbose > 0  && strcmp(pack.saveFigFolder,pack.saveResFolder)
            disp('   Matlab and figure files are in the same folder and will be zipped together.')
        elseif pack.verbose > 1
            disp('   Zipping figure folder...')
        end

        Zip_Delete_Folder(pack.saveFigFolder,(pack.deleteAfterZip == 1 || pack.deleteAfterZip == 3))
    end

    if ~sameFol && (~isempty(pack.saveResFolder) && (pack.zipFiles == 2 || pack.zipFiles == 3))
        if pack.verbose > 1
            disp('   Zipping results folder...')
        end

        Zip_Delete_Folder(pack.saveResFolder,(pack.deleteAfterZip == 2 || pack.deleteAfterZip == 3))
    end

    cd(origPath)
end

end
