function DV_Res = DrainageVolc_Analysis(pack)
%%
% Name: DrainageVolc_Analysis
% Initial Date: 02/04/2021 (mm/dd/yyyy)
% Recent Update Date: 08/20/2025 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to analyze the basin morphometry metrics of a volcanic
% edifice. Including typical topography-based metrics (slope, curvature,
% roughness), drainage-based metrics (drainage area, Hack's Law
% relationship), and channel-based metrics (drainage density, chi). 
%
% NOTE:This script uses TopoToolbox (https://github.com/wschwanghart/topotoolbox), 
%   as well as some other community-created functions (e.g., fitellipse, 
%   ll2utm) and requires Matlab's Mapping, Parallel Computing, and Image 
%   Processing Toolboxes.
%
% This script is compatible with both TopoToolbox 2 and TopoToolbox 3.
%
% ALSO: This script is continuously being updated. Some functionality may  
%   change in future releases. Users are free to use and adapt this script 
%   as needed.
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
%       dx: Grid resolution (in m) of the map.
%
%       hypsIter: Iteration value (0-1) to use for hypsometric/CDF 
%           analysis.
%       roughnessWindows: Window sizes (in m) to calculated roughness. This
%           can be any length, but only four will plot.
%       roughnessType: String describing the type of roughness analysis,
%           values given in TopoToolbox's 'roughness' script.
%       slopeVarianceWindows: Window sizes (in m) to calculate slope
%           variance. This can be any length, but only four will plot.
%       
%       contourSinuosity_ContIter: Iteration value for contour sinuosity
%           analysis. If the value is > 0, it is treated as a topographic
%           contour iteration. If it is a negative value between -1 and 0,
%           it is used as a percentage of the total edifice height.
%       basinContIter: Contour iteration to collect the number of basins
%           per contour length. If given as a negative value between -1 and 
%           0, value will be used as a percentage of the total edifice 
%           height.
%       basinRadIter: Distance interval for number of basins per distance
%           from edifice peak. If given as a negative value between -1 and 
%           0, value will be used as a percentage of the maximum edifice 
%           radius.
%
%       basinTopN: Number of largest basins to analyze for both the
%           slope-area and channel concavity analyses. If basinTopN is a
%           positive integer, the algorithm will analyze that number of
%           largest basins. If basinTopN is between 0 and 1, the algorithm 
%           will analyze basins that fit within the given top percentile of 
%           all basins. If basinTopN is negative and between -1 and 0, the 
%           algorithm will analyze basins that fit within the given top 
%           percentile of elevation.
%       basinTopNFromContLength: Flag to determine if basinTopN elevation
%           percentile should be based on number of basins per contour
%           length value.
%
%       limitHacksLaw: Flag for whether Hack's Law should only be limited
%           to basins with a drainage area greater than the channelization
%           threshold.
%       bsinStatThreshold: Drainage area threshold (in m^2) for whether
%           basin statistics should be conducted.
%       smoothBasinPointWavelength: Wavelength to smooth mid-basin points,
%           to help correct noisy data that can skew basin length values.
%
%       channelThreshold: Drainage area threshold (in m^2) of river 
%           channels.
%       dynamicThreshold: Flag to dynamically calculate the drainage area
%           threshold for channel distinction.
%       dynamicThresholdPixelStep: Pixel step count to calculate the
%           drainage area threshold for channels.
%
%       conformityWavelength: Wavelength to filter topography for use in
%           calculating channel conformity indices (Black et al., 2017). 
%           Setting to NaN will make the wavelenth equivilent to the 
%           edifice's effective radius. 
%       conformityStreamDist: Channel distance over which orientatations
%           and flow paths are calculated for the conformity index.
%       knickpointTolerance: Elevation tolerance value to determine
%           knickpoints using TopoToolbox's knickpointfinder. Lower values
%           determine more knickpoints.
%       MN: M/N value for \Chi analysis. Set as NaN for the script to
%           determine best-fit value.
%       chi_Zcutoff = Cutoff elevation for \Chi analysis. If it is set to 
%           NaN, the elevation is set to the highest boundary point. If set
%           to an imaginary number > 0 the elevation is set to that amount
%           of relief above the lowest edifice elevation (e.g., 50i will
%           set the elevation to 50 m above the lowest point). If set to an
%           imaginary number between -1 and 0, the cutoff is determined as
%           a percent of the total relief above the lowest point.
%       chi_removeUpperBasins = Flag for whether basins that have outlets 
%           above the elevation cutoff should be removed from the \Chi 
%           analysis.
%       concavityType: Type of regression analysis to calculate channel
%           concavity and steepness index. 'ls' is nonlinear least squares
%           from TopoToolbox; 'lad' is absolute deviation from TopoToolbox;
%           'logtrls' is linear least squares from TopoToolbox; 'lin' is
%           linear polygon fit of slope-area logarithmic values; and 'nlin'
%           is nonlinear fit of slope-area values. 'lad' is default.
%           'logtrls' and 'lin' are not recommended.
%
%       Analyze_Divides: Flag for whether divides metrics should be
%           determined.
%       Divide_Order_Cutoff: Threshold value for divide orders to 
%           calculate the Divide Asymmetry Index (DAI).
%       DAI_Integral_BinWidth: Bin width to derive the PDF of DAI values,
%           from which the integral is taken to calculate \Gamma (O'Hara et 
%           al., 2024) 
%
%       parallelProc: Flag for whether parallel processing should be used.
%
%       saveResFolder: Folder to save output structure. Set as empty if not
%           saving results.
%       saveInputs: Flag to save the inputs as a text file.
%
%       plotResults: Flag to plot analysis results.
%       visPlots: Flag to determine whether plots are visible (good for
%           running in background).
%       figPrefix: Prefix to give before saved .fig and .png file names.
%       figTitlePrefix: Prefix for figure titles (e.g, volcano name).
%       saveFigFolder: Folder to save analysis plots (must include '/' or 
%           '\' at the end, and plotResults must be set to 1 to use). Set 
%           to empty if not saving plots.    
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
%   mets: Class variable that contains all analysis results, ordered by 
%   category:
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
%       TopoParams:
%           Slope: Grid of DEM slope values.
%
%           Roughness: Indexed-structure of roughness analysis, includes:
%               Roughness_Grid = Grid of roughness values.
%               WindowSize = Pixel-size of roughness window.
%               ExpectedWindowRes = Expected, user-defined roughness
%                   window.
%               TrueWindowRes = Actual roughness based on grid resolution.
%               Hypsometry_Areas = Normalized area values for roughness.
%               Hypsometry_Values = Roughness CDF values.
%               Hypsometry_NormValues = Nomalized roughness CDF values.
%
%           SlopeVariance_Windows: Indexed-structure of slope variance
%               analysis, includes:
%               SlopeVariance_Grid = Grid of roughness values.
%               WindowSize = Pixel-size of roughness window.
%               ExpectedWindowRes = Expected, user-defined roughness
%                   window.
%               TrueWindowRes = Actual roughness based on grid resolution.
%               Hypsometry_Areas = Normalized area values for roughness.
%               Hypsometry_Values = Roughness CDF values.
%               Hypsometry_NormValues = Nomalized roughness CDF values.
%
%       DrainageParams:
%           Hydrology: Structure of hydrology-based TopoToolbox classes.
%               FD: Flow dirction object (FLOWobj).
%               A: Cumulative drainage area (GRIDobj).
%               D: Flow distance (GRIDobj).
%
%           Drainage_Basins: Structure of drainage basin variables.
%               DB: Separated drainage basins (GRIDobj).
%               DBxy: X-Y coordinates of the drainage basin boundaries, 
%                   separated by NaNs.
%               DBxy_cell: X-Y coordinates of the drainage basin  
%                   boundaries, separated as a cell array.
%               Basin_Contour_ContourP_Count_Length_Area: Array containing 
%                   analyzed contour values (in meters and percent relief  
%                   from main flank), the number of basins in each contour, 
%                   the length of the contour, and the area of the contour.
%               Basin_Contour_ContourP_Count_Length_Area_Channelized: same 
%                   as Basin_Contour_ContourP_Count_Length_Area, but for 
%                   basins with areas greater than the channelization
%                   threshold.
% 
%           Basin_Roughness_Stats: Structure of roughness statistics for
%               each basin. Includes:
%               BasinIDs: Array of basin IDs.
%               Windows: Array of distances over which roughness was
%                   calcualted (i.e., windows), corrected for cell size.
%               Means: Array of roughness means, rows correspond to
%                   basin IDs, columns to windows.
%               Medians: Array of roughness medians, rows correspond to
%                   basin IDs, columns to windows.
%               Stds: Array of roughness standard deviations, rows 
%                   correspond to basin IDs, columns to windows.
%               Mins: Array of roughness minima, rows correspond to
%                   basin IDs, columns to windows.
%               Maxes: Array of roughness maxima, rows correspond to
%                   basin IDs, columns to windows.
%
%           Basin_SlopeVariance_Stats: Structure of slope variance 
%               statistics for each basin. Includes:
%               BasinIDs: Array of basin IDs.
%               Windows: Array of distances over which slope variance was
%                   calcualted (i.e., windows), corrected for cell size.
%               Means: Array of slope variance means, rows correspond to
%                   basin IDs, columns to windows.
%               Medians: Array of slope variance medians, rows correspond 
%                   to basin IDs, columns to windows.
%               Stds: Array of slope variance standard deviations, rows 
%                   correspond to basin IDs, columns to windows.
%               Mins: Array of slope variance minima, rows correspond to
%                   basin IDs, columns to windows.
%               Maxes: Array of slope variance maxima, rows correspond to
%                   basin IDs, columns to windows.
%
%           Basin_Contour_Sinuosity: Structure of contour sinuosity
%               statistics foreach basin. Includes:
%               BasinIDs: Array of basin IDs.
%               Contours: Array of elevations (in m) that define contours.
%               Norm_Contours: Array of normalized elevations (0-1) that
%                   define contours.
%               Individual_Sinuosities: Individual sinuosity values for
%                   each basin and contour. Rows correspond to basins,
%                   columns to elevations.
%               Means: Array of average sinuosity values for each basin.
%               Medians: Array of median sinuosity values for each basin.
%               Stds: Array of sinuosity value standard deviations for each
%                   basin.
%               Mins: Array of minimum sinuosity values for each basin.
%               Maxes: Arry of maximum sinuosity values for each basin.
%
%           Top_Drainage_Basins: Structure of  highest-ranking basins
%               Basin_TopN: Value used to define the top basins.
%               TopN_basinIDs: Basin IDs of the top basins.
%               TopDBxy: X-Y coordinates of the largest drainage basin 
%                   boundaries, separated by NaNs.
%               TopDBxy_cell: X-Y coordinates of the largest drainage basin 
%                   boundaries, separated as a cell array.
%               TopN_All_Area_Slope: Drainage area and slope values of the
%                   largest N basins. Used for slope-area analysis.
%               TopN_Flow_Area_Slope: Drainage area and slope values of 
%                   only longest flow paths of the top N basins. Used for 
%                   slope-area analysis.  
%               TopN_AreaThreshold_Rs: Array that contains the
%                   dynamically-calculated drainage area threshold for
%                   channels, the r^2 value of the hillslope regression, 
%                   r^2 value of channel regression, and r' value used to 
%                   determine the threshold.
%               TopN_TransitionThreshold_A: Area threhold for transition 
%                   zone initiation.
%               TopN_MN: M/N value derived from linear regression of 
%                   slope-area  plots to determine best drainage area 
%                   threshold for channelization.
%
%           Hacks_Law: Structure of Hack's Law analysis values
%               HackLawFit_BasinLength: Best-fit parameters for Hack's Law 
%                   (relationship between drainage area and length) based  
%                   on basin length.
%               HackLawFit_FlowLength: Best-fit parameters for Hack's Law 
%                   based on basin flow length.
%               HackLawDeviation_BasinLength: Array of basin deviations  
%                   from best-fitting Hack's Law parameters using basin 
%                   length. Columns are basin ID, basin area, basin length, 
%                   basin length deviation from best-fit Hack parameters, 
%                   and flag for whether the basin was used to derive the 
%                   best-fit power law.
%               HackLawDeviation_FlowLength: Array of basin deviations from 
%                   best-fitting Hack's Law parameters using flow length.
%                   Columns are basin ID, basin area, basin length, basin
%                   length deviation from best-fit Hack parameters, and 
%                   flag for whether the basin was used to derive the 
%                   best-fit power law.
%               HackLawDeviation_BasinLength_Map: GRIDobj of basin  
%                   deviations from Hack's Law derived from basin lengths.
%               HackLawDeviation_FlowLength_Map: GRIDobj of basin  
%                   deviations from Hack's Law derived from flow lengths.
%               HacksLawDeviation_Titles: Cell array of titles for Hack's
%                   Law Deviation arrays.
%
%           Basin_Statistics: Structure of basin statistic variables.
%               Basin_Statistics: Matrix of basin IDs, drainage areas, max 
%                   flow lengths, basin lengths, widths, reliefs, 
%                   orientations, hypsometry integrals, mean slopes, 
%                   sinuosity, Euclidean distance between headwater and 
%                   outlet, along-basin distance from headwater to largest 
%                   width, and slope variance.
%               Basin_Statistics_Grids: Structure of grids for overall 
%                   basin lengths, heights, widths, orientations, flow 
%                   lengths, hypsometry integrals, mean slopes, s
%                   inuosities, drainage  area, and slope variance.
%               Basin_Cross_Statistics: Array containing cross-basin 
%                   stistics. Rows are ordered as: basin ids, drainage  
%                   sampling x-coordinates, drainage sampling y-coordinates, 
%                   drainage sampling elevations, cross-basin widths 
%                   (perpendicular to broad basin orientation) cross-basin 
%                   relief (perpendicular to broad basin orientation), and 
%                   incision ratio  (relief / width).
%               Basin_Statistics_Titles: Cell array of metrics titles
%                   corresponding with Basin_Statistics columns.
%
%           Radial_Analysis: Structure of basins along radial distance from
%                   the edifice's peak.
%               Basin_Count: Matrix of basin count results. Columns are bin 
%                   distance, normalized bin distance, all basin count, and 
%                   basin count of only those above areaThreshold.
%               Basin_Count_Titles: Cell array of table columns.
%               Radial_Distances: GRIDobj of radial distances from the 
%                   edifice's peak.
%               Normalized_Radial_Distances: GRIDobj of raidal distances 
%                   from the edifice's peak, normalized by the maximum 
%                   radius.
%               Basin_IDs: Basin ID's at each interval.
%               Basin_IDs_Channelized: Basin ID's above the channelization 
%                   threshold at each interval.
%
%           Hypsometry: Structure of basin hypsometry calculations.
%               DB_Hyps_Topo: Normalized topography values for drainage  
%                   basin CDF.
%               DB_Hyps_numDB: Drainage basin CDF values.
%               DB_Hyps_numDB_Channelized: CDF values of drainage basins
%                   greater than the channelization threshold.
%
%       ChannelParams:
%           Channels: Structure of variables describing identified channels.
%               ChannelThreshold: Drainage area threshold used for channel
%                   distinction.
%               S: Determined river channels (STREAMobj).
%               Stream_Order: Array of Strahler stream orders for each 
%                   segment of S (needs to be paired with S).
%               Max_Stream_Order: GRIDobj of maximum stream orders in each
%                   basin.
%               Knickpoint_ID_BID_XY_StreamDist_DZ: Array of knickpoints
%                   determined from TopoToolbox's knickpoint finder. 
%                   Columns are unique knickpoint ID, associated basin ID, 
%                   knickpoint XY coordinate, upstream distance to 
%                   knickpoint, and the elevation difference at the 
%                   knickpoint (i.e., the knickpoint's magnitude).
%
%           Drainage_Density: Structure of drainage density variables.
%               DD: Channel drainage densities (needs to be paired with S).
%               MaxDD: Maximum drainage density of each basin (GRIDobj).
%               BasinDD: Total drainage density of each basin (GRIDobj).
%               TotalDD: Cumulative drainage density of entire edifice.
%
%           Conformity: Structure of channel conformity. The conformity
%                   index (Black et al., 2017) calcuates the azimuthal
%                   difference between a channel segment and the flow path
%                   derived from filtered topography, thus analyzing how 
%                   well a channel conforms to longer wavelength 
%                   topography. Afterwards, the cosine of the difference is 
%                   calculated to put the value between -1 and 1. Here, 
%                   only absolute values are considered. 
%               Mean_Total_Conformity: Mean conformity value for every
%                   channel segment of the edifice.
%               Mean_Basin_Conformity: Array of mean conformity values for
%                   each basin. First column is basin ID, second column is
%                   mean conformity.
%               Mean_Basin_Conformity_Map: GRIDobj of mean conformity
%                   values for each basin.
%               Mean_Basin_StrahlerOrder_Conformity: Structured array of
%                   mean conformity values broken down by both basin and 
%                   Strahler order. Fields include:
%                   Basin_ID: Basin IDs;
%                   Stream_Order: Array of stream orders in the basin.
%                   Mean_Conformities: Array of mean conformity value for
%                       each Strahler order.
%               Segment_Information: Structured array that contains all the
%                   information for how conformity was calculated for each 
%                   segment. Includes:
%                   Basin_ID: Basin ID.
%                   X: Channel segment x-coordinates.
%                   Y: Channel segment y-coordinates.
%                   Z: Channel segment elevations.
%                   Dist: Channel distances along the segment.
%                   Order: Strahler order of the segment.
%                   Orientation: Segment azimuthal orientation.
%                   LongWavelength_EndPoints: End points of a flow path the
%                       same length as the channel segment calculated over
%                       the filtered topography.
%                   LongWavelength_Orientation: Azimuthal orientation of
%                       the fitered topography flowpath.
%                   Orientation_Difference: Absolute difference between
%                       channel segment and flow path orientations.
%                   Conformity: Cosine of Orientation_Difference.
%               Filtered_Topography: GRIDobj of the filtered topography.
%
%           Channel_Concavity: Structure of channel concavity variables.
%               Concavity_Streams: Streams used to analyze channel 
%                   concavity.
%               Concavity_Stats: Stream concavity statistics, including
%                   contributing drainage area, gradient, ks, and best-fit
%                   concavity index (theta)
%               Concavity_BasinIDs: Basin IDs used for analysis.
%               Concavity_Mean: Mean concavity index of all analyzed 
%                   streams.
%               Concavity_DEM: DEM of channel concavities (GRIDobj).
%
%           Chi: Structure of channel chi values.
%               Total_Chi: Structure of channel Chi values, calculated by
%                       assuming all basins have the same M/N value (either  
%                       given by the user, or determined by TopoToolbox). 
%                   Chi_S: STREAMobj used to calculate Chi.
%                   Chi_S_TopN: STREAMobj of only basins that meet the
%                       basinTopN criteria. Used to calculate M/N if it is 
%                       not supplied.
%                   MN: M/N value either supplied by the user or determined
%                       from Chi_S_TopN.
%                   Chi: Array of Chi values, to be paired with Chi_S.
%                   Maximum_Chi: GRIDobj of maximum Chi value for each 
%                       basin.
%                   Upstream_Chi: GRIDobj of Chi values projected upslope 
%                       from the channel to divides. Used to compare Chi  
%                       across divides and determine divide mobility. 
%               Basin_Chi: Structure of channel Chi values, calculated by
%                       assuming each basin has its own M/N value, as 
%                       determined by TopoToolbox).
%                   Chi_S: Structured array of individual basin STREAMobjs.
%                           Contains:
%                       BasinID: Basin ID.
%                       S: Basin STREAMobj.
%                   Chi: Structured array of individual basin Chi values. 
%                           Contains: 
%                       BasinID: Basin ID:
%                       MN: Best-fitting M/N value.
%                       Chi: Array of Chi values, to be paired with Chi_S.
%                   Maximum_Chi: GRIDobj of maximum Chi value for each 
%                       basin.
%                   Upstream_Chi: GRIDobj of Chi values projected upslope 
%                       from the channel to divides. Used to compare Chi  
%                       across divides and determine divide mobility. 
%               Chi_Cutoff_Elevation: Cutoff elevtion used to calculate Chi.
%
%       DivideParams:
%           DivideTopo: Divide catalog ordered by topography.
%               DVD: Divide structure (DIVIDEobj), with ordering.
%               X: Divide X location.
%               Y: Divide Y location.
%           Divide_ReliefDistance: List of distances to common streams across
%               divides.
%           VerticalDistance: Grid of vertical distances between divides
%               and channels (GRIDobj).
%           Divide_AsymmetryIndex: Asymmetry index from Scherler & 
%               Schwanghart (2020)
%           DAI_Gamma: Structure of values used to calculate Gamma. Gamma
%               is the integral of the PDF of Divide Asymmetry Index
%               values (O'Hara et al., 2024). Contains:
%               Bins: DAI bins used for the PDF.
%               Bin_Midpoints: Midpoint values of the bins.
%               PDF: PDF values.
%               Gamma: Integral of the PDF.
%           DAI_Gamma_Corrected: Structure of values to calcualte Gamma
%               with with DAI values of 1 removed to create more consistent
%               values. Has the same structure as DAI_Gamma.
%           Divide_TotalChiDifference: Differences in Chi across divides, 
%               using the same M/N value for all basins.
%           Divide_BasinChiDifference: Differences in Chi across divides, 
%               using different, best-fit M/N values for each basin.
%           Junction_X_Y_C_Z_D_A: Divide junction information, including
%               X-coordinates, Y-coordinates, conjuction index, elevation,
%               distance along divide, and asymmetry index.

%% Unpack input
pack = DrainageVolc_DefaultVals(pack);
if pack.verbose > 0
    disp(sprintf('USING DRAINAGEVOLC VERSION %s',pack.version))
    disp('Unpacking Input and Importing Data...')
end

dx = pack.dx;
basinTopN = pack.basinTopN;
channelThreshold = pack.channelThreshold;
MN = pack.MN;

DV_Res.GeneralParams.Version = pack.version;
DV_Res.GeneralParams.StartTime = datetime('now');

DrainageVolc_CheckIters(pack.hypsIter,pack.basinContIter,pack.basinRadIter,basinTopN)

%% Import Shapefiles
[boundaryXY,craterXY,maskXY] = Import_Shapefiles(pack.boundaryXY,pack.craterXY,pack.maskXY,pack.verbose);

DV_Res.GeneralParams.inputs = pack;
DV_Res.GeographicParams.boundaryXY = boundaryXY;
DV_Res.GeographicParams.craterXY = craterXY;
DV_Res.GeographicParams.maskXY = maskXY;

%% Load and Cut DEM
[DEM0,~,~,~,DEM,X,Y,Z] = Import_DEM(pack.tifFile,dx,1,boundaryXY,craterXY,maskXY,pack.verbose);
[~,~] = meshgrid(X,Y);

DV_Res.GeographicParams.DEM0 = DEM0;
DV_Res.GeographicParams.DEM = DEM;

%% Basic Topography Metrics
if pack.verbose > 0
    disp('Calculating Topography Metrics...')
    disp('   Slope...')
end

slp = gradient8(DEM,'degree');
DV_Res.TopoParams.Slope = slp;

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

    slVarGrid = CalculateSlopeVarianceWindow(slp,slVarI);

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

DV_Res.TopoParams.SlopeVariance_Windows = slopeVarianceVals;

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

DV_Res.TopoParams.Roughness = roughnessVals;

%% Drainage Metrics
if pack.verbose > 0
    disp('Calculating Drainage Metrics...')
end
FD = FLOWobj(DEM,'preprocess','none');

DV_Res.DrainageParams.Hydrology.FD = FD;

% Drainage Area
if pack.verbose > 0
    disp('   Drainage Area...')
end
A = flowacc(FD);
A.Z(isnan(DEM.Z)) = NaN;

% Flow Distance
if pack.verbose > 0
    disp('   Flow Distance...')
end
D = flowdistance(FD);
D.Z(isnan(DEM.Z)) = NaN;

% Drainage Basin
if pack.verbose > 0
    disp('   Drainage Basins...')
end
DB = drainagebasins(FD);
DB = shufflelabel(DB);
DB.Z = single(DB.Z);
DB.Z(isnan(DEM.Z)) = NaN;

[DBxy,DBxy_cell] = DrainageVolc_GetBasinBoundaries(DB,A,dx,channelThreshold);

DV_Res.DrainageParams.Hydrology.A = A;
DV_Res.DrainageParams.Hydrology.D = D;

DV_Res.DrainageParams.Drainage_Basins.DB = DB;
DV_Res.DrainageParams.Drainage_Basins.DBxy = DBxy;
DV_Res.DrainageParams.Drainage_Basins.DBxy_cell = DBxy_cell;

% Drainage Basin Stats
if pack.verbose > 0
    disp('   Drainage Basin Statistics...')
end

[tot_db_stats,tot_db_grids,cross_db_stats] = DrainageVolc_Collect_DrainageBasin_Stats_wCross(Z,slp,A,DB,FD,D,pack.parallelProc,pack.smoothBasinPointWavelength,pack.basinStatThreshold,pack.verbose);

% Drainage Basins Per Contour
Basin_Contour_ContourP_Count_Length_Area = DrainageVolc_Determine_Basin_Per_Contour(DEM,DB,0,pack.basinContIter,NaN);
DV_Res.DrainageParams.Drainage_Basins.Basin_Contour_ContourP_Count_Length_Area = Basin_Contour_ContourP_Count_Length_Area;

% Roughness and Windowed Slope Variance
if pack.verbose > 0
    disp('   Drainage Basin Winowed Roughness and Slope Variance Statistics...')
end
[db_roughnessVals,db_svVals] = DrainageVolc_Collect_DB_Roughness_SlpVar(DB,roughnessVals,slopeVarianceVals,pack.verbose);
DV_Res.DrainageParams.Basin_Roughness_Stats = db_roughnessVals;
DV_Res.DrainageParams.Basin_SlopeVariance_Stats = db_svVals;

% Contour Sinuosity
if pack.verbose > 0
    disp('   Drainage Basin Contour Sinuosity Values...')
end
 db_contSin = DrainageVolc_Collect_DB_ContourSinuosity(DEM,DB,pack.contourSinuosity_ContIter,pack.verbose);
DV_Res.DrainageParams.Basin_Contour_Sinuosity_Stats = db_contSin;

% Drainage Area - Slope Plots
if pack.verbose > 0
    disp('   Slope-Area Relationship...')
end
if pack.basinTopNFromContLength
    tmp = Basin_Contour_ContourP_Count_Length_Area;
    tmp(tmp(:,2) >= .9,:) = [];
    contLength = tmp(:,3)./tmp(:,4);
    diffcontLength = [0;diff(contLength)];

    ii = find(diffcontLength<0,1,'last');
    basinTopN = (tmp(ii,2))-1;
end
[topNAS,topNFlowAS,topDBxy,DBxy_cell,topN_basinIDs] = DrainageVolc_Collect_SlopeAreas(tot_db_stats,DB,D,FD,slp,A,dx,basinTopN,DEM);
DV_Res.DrainageParams.Top_Drainage_Basins.TopDBxy = topDBxy;
DV_Res.DrainageParams.Top_Drainage_Basins.TopDBxy_cell = DBxy_cell;
DV_Res.DrainageParams.Top_Drainage_Basins.TopN_All_Area_Slope = topNAS;
DV_Res.DrainageParams.Top_Drainage_Basins.TopN_Flow_Area_Slope = topNFlowAS;
DV_Res.DrainageParams.Top_Drainage_Basins.Basin_TopN = basinTopN;
DV_Res.DrainageParams.Top_Drainage_Basins.TopN_basinIDs = sortrows(topN_basinIDs);

% Find Dynamic Threshold
if pack.dynamicThreshold
    [basinAreaThreshold,transitionThreshold,MN] = DrainageVolc_SA_PiecewiseRegression(topNFlowAS,pack.dynamicThresholdPixelStep*dx^2);
    channelThreshold = basinAreaThreshold(1);
    
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_AreaThreshold_Rs = basinAreaThreshold;
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_TransitionThreshold_A = transitionThreshold;
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_MN = abs(MN);
else
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_AreaThreshold_Rs = [];
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_TransitionThreshold_A = [];
    DV_Res.DrainageParams.Top_Drainage_Basins.TopN_MN = [];
end
DV_Res.ChannelParams.Channels.ChannelThreshold = channelThreshold;

% Hack's Relationship (L = c*A^b) Using Basin Length
if pack.verbose > 0
    disp('   Hack''s Law Parameterization...')
end

if pack.limitHacksLaw
    [tmpPF,PF,lDev_Map,LDev_Titles] = DrainageVolc_HacksLaw(DB,tot_db_stats(:,1),...
        tot_db_stats(:,2),tot_db_stats(:,4),channelThreshold);
else
    [tmpPF,PF,lDev_Map,LDev_Titles] = DrainageVolc_HacksLaw(DB,tot_db_stats(:,1),...
        tot_db_stats(:,2),tot_db_stats(:,4),NaN);
end

DV_Res.DrainageParams.Hacks_Law.HackLawFit_BasinLength = PF;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_BasinLength = tmpPF;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_BasinLength_Map = lDev_Map;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_Titles = LDev_Titles;

% Hack's Relationship (L = c*A^b) Using Flow Length
if pack.limitHacksLaw
    [tmpPF,PF,lDev_Map,LDev_Titles] = DrainageVolc_HacksLaw(DB,tot_db_stats(:,1),...
        tot_db_stats(:,2),tot_db_stats(:,3),channelThreshold);
else
    [tmpPF,PF,lDev_Map,LDev_Titles] = DrainageVolc_HacksLaw(DB,tot_db_stats(:,1),...
        tot_db_stats(:,2),tot_db_stats(:,3),NaN);
end
DV_Res.DrainageParams.Hacks_Law.HackLawFit_FlowLength = PF;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_FlowLength = tmpPF;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_FlowLength_Map = lDev_Map;
DV_Res.DrainageParams.Hacks_Law.HackLawDeviation_Titles = LDev_Titles;

DV_Res.DrainageParams.Basin_Statistics.Basin_Statistics = tot_db_stats;
DV_Res.DrainageParams.Basin_Statistics.Basin_Statistics_Grids = tot_db_grids;
DV_Res.DrainageParams.Basin_Statistics.Basin_Cross_Statistics = cross_db_stats;
DV_Res.DrainageParams.Basin_Statistics.Basin_Statistics_Titles = {'ID','Drainage Area','Flowpath Length',...
    'Basin Length','Basin Width','Basin Height','Basin Orientation',...
    'Basin Hypsometry Integral','Mean Basin Slope','Flowpath Sinuosity','Basin Euclidean Distance',...
    'Basin Length Distance to Largest Width','Slope Variance'};

% Drainage Basins Per Radial Distance
if pack.verbose > 0
    disp('   Drainage Basin by Radial Distance...')
end
[BCount_Table,BCount_Titles,BIDs,BIDs_AboveA,R_DEM,RNorm_DEM] = DrainageVolc_CalculateRadialDistanceNumbers(DEM,DB,tot_db_grids.DrainageArea,pack.basinRadIter,channelThreshold);
DV_Res.DrainageParams.Radial_Analysis.Basin_Count = BCount_Table;
DV_Res.DrainageParams.Radial_Analysis.Basin_Count_Titles = BCount_Titles;
DV_Res.DrainageParams.Radial_Analysis.Radial_Distances = R_DEM;
DV_Res.DrainageParams.Radial_Analysis.Normalized_Radial_Distances = RNorm_DEM;
DV_Res.DrainageParams.Radial_Analysis.Basin_IDs = BIDs;
DV_Res.DrainageParams.Radial_Analysis.Basin_IDs_Channelized = BIDs_AboveA;

% Drainage Basin Hypsometry Counts
[DB_Hyps_Topo,DB_Hyps_numDB,DB_Hyps_numDB_aboveA] = DraingeVolc_Count_VolcTopoBasins(DEM.Z,DB.Z,tot_db_grids.DrainageArea.Z,pack.hypsIter,channelThreshold);
DV_Res.DrainageParams.Hypsometry.DB_Hyps_Topo = DB_Hyps_Topo;
DV_Res.DrainageParams.Hypsometry.DB_Hyps_numDB = DB_Hyps_numDB;
DV_Res.DrainageParams.Hypsometry.DB_Hyps_numDB_Channelized = DB_Hyps_numDB_aboveA;

% Channelized Drainage Basins Per Contour
Basin_Contour_ContourP_Count_Length_Area_AboveA = DrainageVolc_Determine_Basin_Per_Contour(DEM,DB,tot_db_grids.DrainageArea,pack.basinContIter,channelThreshold);
DV_Res.DrainageParams.Drainage_Basins.Basin_Contour_ContourP_Count_Length_Area_Channelized = Basin_Contour_ContourP_Count_Length_Area_AboveA;

%% Channel Metrics
if pack.verbose > 0
    disp('Calculating Channel Metrics...')
end
S = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');
DV_Res.ChannelParams.Channels.S = S;

if isempty(S)
    warning('NO CHANNELS FOUND, DECREASE DRAINAGE AREA THRESHOLD')
end

if pack.verbose > 0
    disp('   Strahler Order...')
end

SO = streamorder(S);
DV_Res.ChannelParams.Channels.Stream_Order = SO;

tmpSOG = STREAMobj2GRIDobj(S,SO);
[tmpSOG,~,~] = GRIDobj2mat(tmpSOG);

% Drainage density
if pack.verbose > 0
    disp('   Drainage Density...')
end
DD = drainagedensity(S,FD,'nal');
DDG = drainagedensity(S,FD,'grid');
DB_DD = DB;
DB_DD.Z(:) = NaN;
DB_mDD = DB_DD;
uniDB = unique(DB.Z(:));
uniDB(isnan(uniDB)) = [];

sTmp = STREAMobj2GRIDobj(S);
SOG = DB;
SOG.Z(:) = NaN;

for i = 1:length(uniDB)
    tmpDB = DB.Z == uniDB(i);
    maxVal = max(max(DDG.Z(tmpDB==1)));
    DB_mDD.Z(tmpDB==1) = maxVal;

    maxSO = max(max(tmpSOG(tmpDB==1)));

    tmpZ = DEM;
    tmpZ.Z(DB.Z~=uniDB(i)) = NaN;
    tmpZ.Z(~sTmp.Z) = NaN;
    [ii,jj] = find(tmpZ.Z==min(tmpZ.Z(:)),1);
    if ~isempty(ii)
        DB_DD.Z(tmpDB==1) = DDG.Z(ii,jj);
        SOG.Z(tmpDB==1) = maxSO;
    end
end

DV_Res.ChannelParams.Channels.Max_Stream_Order = SOG;

DV_Res.ChannelParams.Drainage_Density.DD = DD;
DV_Res.ChannelParams.Drainage_Density.MaxDD = DB_mDD;
DV_Res.ChannelParams.Drainage_Density.BasinDD = DB_DD;

DV_Res.ChannelParams.Drainage_Density.TotalDD = DrainageVolc_Collect_Total_DrainageDensity(Z,S,A,DB,dx);

% Conformity Index
if pack.verbose > 0
    disp('   Channel Conformity Index...')
end

[totalConformity,basinConformity,basinConformityMap,basinOrderConformity,segmentInfo,filteredDEM] = DrainageVolc_Calculate_ConfomityIndex(DEM0,DEM,boundaryXY,DB,channelThreshold,pack.conformityWavelength,pack.conformityStreamDist,pack.verbose);
DV_Res.ChannelParams.Conformity.Mean_Total_Conformity = totalConformity;
DV_Res.ChannelParams.Conformity.Mean_Basin_Conformity = basinConformity;
DV_Res.ChannelParams.Conformity.Mean_Basin_Conformity_Map = basinConformityMap;
DV_Res.ChannelParams.Conformity.Mean_Basin_StrahlerOrder_Conformity = basinOrderConformity;
DV_Res.ChannelParams.Conformity.Segment_Information = segmentInfo;
DV_Res.ChannelParams.Conformity.Filtered_Topography = filteredDEM;

% Concavity Analysis
if pack.verbose > 0
    disp('   Channel Concavity...')
end
[topN_Ss,topN_S_Stats] = DrainageVolc_Collect_Stream_Concavity(tot_db_stats,DEM,A,DB,FD,topN_basinIDs,channelThreshold,pack.concavityType);

DV_Res.ChannelParams.Channel_Concavity.Concavity_Streams = topN_Ss;
DV_Res.ChannelParams.Channel_Concavity.Concavity_Stats = topN_S_Stats;
DV_Res.ChannelParams.Channel_Concavity.Concavity_BasinIDs = topN_basinIDs;

conDEM = DEM;
conDEM.Z(:,:) = NaN;
for i = 1:length(topN_basinIDs)
    conDEM.Z(DB.Z==topN_basinIDs(i)) = abs(topN_S_Stats(i).theta);
end

allTheta = [];
for i = 1:size(topN_Ss,1)
    allTheta = [allTheta;topN_S_Stats(i).theta];
end
DV_Res.ChannelParams.Channel_Concavity.Concavity_Mean = nanmean(allTheta);
DV_Res.ChannelParams.Channel_Concavity.Concavity_DEM = conDEM;

% Chi
%   This function uses Beysian optimization to find the best-fitting m/n
%   ratio for chi.
if pack.verbose > 0
    disp('   Chi Analysis...')
end

[Total_Chi_Class,Individual_Chi_Class,Chi_Cutoff_Elevation] = DrainageVolc_Calculate_Chi(DEM0,DEM,DB,...
    boundaryXY,craterXY,channelThreshold,topN_basinIDs,MN,pack.plotResults,pack.saveFigFolder,pack.figPrefix,pack.parallelProc,pack.chi_Zcutoff,pack.chi_removeUpperBasins,pack.verbose);

DV_Res.ChannelParams.Chi.Total_Chi = Total_Chi_Class;
DV_Res.ChannelParams.Chi.Basin_Chi = Individual_Chi_Class;
DV_Res.ChannelParams.Chi.Chi_Cutoff_Elevation = Chi_Cutoff_Elevation;

% Knickpoint Analysis
if pack.verbose > 0
    disp('   Knickpoint Locations...')
end
kn_ID_BID_XY_Dist_DZ = DrainageVolc_Find_Knickpoints(DEM,DB,channelThreshold,pack.knickpointTolerance,pack.parallelProc,pack.verbose);
DV_Res.ChannelParams.Channels.Knickpoint_ID_BID_XY_StreamDist_DZ = kn_ID_BID_XY_Dist_DZ;

%% Divide Metrics
if ~pack.Analyze_Divides
    DV_Res.DivideParams = [];
else
    if pack.verbose > 0
        disp('Calculating Divide Metrics...')
    end
    % Collect divides
    if isempty(S)
        warning('IGNORING DIVIDES SINCE NO CHANNELS FOUND')
        DV_Res.DivideParams.Divide_Topo = [];
        DV_Res.DivideParams.Divide_ReliefDistance = [];
        DV_Res.DivideParams.Divide_AsymmetryIndex = [];
        DV_Res.DivideParams.DAI_Gamma = [];
        DV_Res.DivideParams.DAI_Gamma_Corrected = [];
        DV_Res.DivideParams.Divide_TotalChiDifference = [];
        DV_Res.DivideParams.Divide_BasinChiDifference = [];
        DV_Res.DivideParams.Junction_X_Y_C_Z_D_A = [];
        DV_Res.DivideParams.VerticalDistance = [];
    else
        DVD_Topo = DIVIDEobj(FD,S,'type','topo');
        DVD_Topo = cleanedges(DVD_Topo,FD);
        DVD_Topo = sort(DVD_Topo);
        DVD_Topo = divorder(DVD_Topo);
        
        % Order divides and find their roots
        [X_Topo,Y_Topo] = ind2coord(DVD_Topo,DVD_Topo.IX(1:end-1));
        
        Divide_Topo.DVD = DVD_Topo;
        Divide_Topo.X = X_Topo;
        Divide_Topo.Y = Y_Topo;
        DV_Res.DivideParams.Divide_Topo = Divide_Topo;
        
        % Calculate Divide Asymmetry Index
        if pack.verbose > 0
            disp('   Divide Asymmetry Index (DAI)...')
        end
        vertDist = vertdistance2stream(FD,S,DEM);
        vertDist.Z(isinf(vertDist.Z)) = NaN;
        
        [pp,qq] = getvalue(DVD_Topo,vertDist);
        reliefDistance = abs(diff(abs([pp,qq]),1,2));
        asymIndex = reliefDistance./sum([pp,qq],2);
        asymIndex(asymIndex<0) = 0;
        asymIndex(isinf(asymIndex)) = NaN;
        
        DV_Res.DivideParams.Divide_ReliefDistance = reliefDistance;
        DV_Res.DivideParams.Divide_AsymmetryIndex = asymIndex;
    
        % Calculate DAI integral
        if pack.verbose > 0
            disp('   DAI Integral (Gamma)...')
        end
    
        tmp = asymIndex;
        tmp(DVD_Topo.order < pack.Divide_Order_Cutoff) = NaN;
        tmp(isnan(tmp)) = [];
    
        asBins = 0:pack.DAI_Integral_BinWidth:1;
        [DAI_PDF,~] = histcounts(tmp,asBins,'Normalization','pdf');
        gamma = trapz(pack.DAI_Integral_BinWidth,DAI_PDF);
        gammaStruct.PDF = DAI_PDF;
        gammaStruct.Bins = asBins;
        gammaStruct.Bin_Midpoints = asBins(1:end-1)+pack.DAI_Integral_BinWidth/2;
        gammaStruct.Gamma = gamma;
    
        DV_Res.DivideParams.DAI_Gamma = gammaStruct;
    
        % Corrected DAI integral (with no 1s)
        tmp(tmp==1) = [];
        asBins = 0:pack.DAI_Integral_BinWidth:1;
        [DAI_PDF,~] = histcounts(tmp,asBins,'Normalization','pdf');
        gamma = trapz(pack.DAI_Integral_BinWidth,DAI_PDF);
        gammaStruct2.PDF = DAI_PDF;
        gammaStruct2.Bins = asBins;
        gammaStruct2.Bin_Midpoints = asBins(1:end-1)+pack.DAI_Integral_BinWidth/2;
        gammaStruct2.Gamma = gamma;
    
        DV_Res.DivideParams.DAI_Gamma_Corrected = gammaStruct2;
    
        % Calculate Chi differences across the divide
        if pack.verbose > 0
            disp('   Chi differences...')
        end
        try
            [pp,qq] = getvalue(DVD_Topo,Total_Chi_Class.Upstream_Chi);
            totalChiDifference = abs(diff(abs([pp,qq]),1,2));
        catch
            totalChiDifference = NaN;
        end
        DV_Res.DivideParams.Divide_TotalChiDifference = totalChiDifference;
    
        try
            [pp,qq] = getvalue(DVD_Topo,Individual_Chi_Class.Upstream_Chi);
            basinChiDifference = abs(diff(abs([pp,qq]),1,2));
        catch
            basinChiDifference = NaN;
        end
        DV_Res.DivideParams.Divide_BasinChiDifference = basinChiDifference;
        
        % Determine Divide Connectivity
        if pack.verbose > 0    
            disp('   Divide Connectivity...')
        end
        try
            [CJ,CJ_x,CJ_y] = jctcon(DVD_Topo);
            Junction_X_Y_C_Z_D_A = DrainageVolc_Collect_JunctionStats(DEM,DVD_Topo,asymIndex,CJ,CJ_x,CJ_y);
        catch er
            Junction_X_Y_C_Z_D_A = [NaN,NaN,NaN,NaN,NaN,NaN];
        end
        
        DV_Res.DivideParams.Junction_X_Y_C_Z_D_A = Junction_X_Y_C_Z_D_A;
        DV_Res.DivideParams.VerticalDistance = vertDist;
    end
end

%% Save All Results
DV_Res.GeneralParams.EndTime = datetime('now');
if ~isempty(pack.saveResFolder)
    if pack.verbose > 0   
        disp('Saving Results...')
    end
    save([pack.saveResFolder,pack.figPrefix,'DrainageVolc_Results.mat'],'DV_Res')
end

%% Save Inputs
if pack.saveInputs && ~isempty(pack.saveResFolder)
    if pack.verbose > 0   
        disp('Writing Input Text File...')
    end

    tmpTif = pack.tifFile;
    tmpMask = pack.maskXY;
    tmpCrat = pack.craterXY;

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

    Export_Inputs(pack,'DrainageVolc');

    pack.tifFile = tmpTif;
    pack.maskXY = tmpMask;
    pack.craterXY = tmpCrat;
end

%% Plot Results
if pack.plotResults
    if pack.verbose > 0   
        disp('Plotting Results...')
    end
    DrainageVolc_Plots(DV_Res);
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

