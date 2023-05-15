function mets = DrainageVolc_Analysis(pack)
%%
% Name: DrainageVolc_Analysis
% Date: 02/04/2021 (mm/dd/yyyy)
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
% ALSO: This script is continuously being updated. Some functionality is 
%   implemented in DrainageVolc that needs fixed (e.g, curvature analyses); 
%   this is underway. Users are free to use and adapt this script as needed.
%
% Input: 
%   pack: Strucure of input parameters, includes:
%       tifFile: .tif file to analyze. 
%       boundaryXY: 2xN matrix or shapefile name of edifice boundary 
%           locations. If matrix, first column is X and second column is Y.
%       maskMap: Array, or .shp file and path, of mask regions X- and 
%           Y-coordinates to ignore in analysis.
%       craterXY: 2xN matrix or shapefile name of crater boundary location
%           that is ignored during DEM filling. If matrix, first column is
%           X and second column is Y.
%
%       dx: Grid resolution (in m) of the map.
%       
%       hypsIter: Iteration value (0-1) to use for hypsometric/CDF 
%           analysis.
%       roughnessWindows: Window sizes (in m) to calculated roughness.
%
%       channelThreshold: Drainage area threshold (in m^2) of river 
%           channels.
%       dynamicThreshold: Flag to dynamically calculate the drainage area
%           threshold for channel distinction.
%       dynamicThresholdPixelStep: Pixel step count to calculate the
%           drainage area threshold for channels.
%
%       basinContIter: Contour iteration to collect the number of basins
%           per contour length. If given as a negative value between -1 and 
%           0, value will be used as a percentage of the total edifice 
%           height.
%       basinTopN: Number of largest basins to analyze for both the
%           slope-area and channel concavity analyses. If basinTopN is 
%           between 0 and 1, the algorithm will analyze basins that fit 
%           within the given top percentile of all basins. If basinTopN is
%           negative and between -1 and 0, the algorithm will analyze
%           basins that fit within the given top percentile of elevation.
%       basinTopNFromContLength: Flag to determine if basinTopN elevation
%           percentile should be based on number of basins per contour
%           length value.
%       Analyze_Divides: Flag for whether divides metrics should be
%           determined.
%
%       smoothBasinPointWavelength: Wavelength to smooth mid-basin points,
%           to help correct noisy data that can skew basin length values.
%       parallelProc: Flag for whether parallel processing should be used.
%       MN: M/N value for \Chi analysis. Set as NaN for the script to
%           determine best-fit value.
%
%       saveResFolder: Folder to save output structure. Set as empty if not
%           saving results.
%
%       plotResults: Flag to plot analysis results.
%       figPrefix: Prefix to give before saved .fig and .png file names.
%       saveFigFolder: Folder to save analysis plots (must include '/' or 
%           '\' at the end, and plotResults must be set to 1 to use). Set 
%           to empty if not saving plots.                                         
%
% Output: 
%   mets: Class variable that contains all analysis results, ordered by 
%   category:
%       GeneralParams:
%           DEM0: Initial DEM (GRIDobj).
%           X: Vector of DEM X coordinates.
%           Y: Vector of DEM Y coordinates.
%           Z: Grid of elevations.
%           DEM: Analyzed DEM (GRIDobj).
%           boundaryXY: Volcano boundary X- and Y-coordinates.
%           craterXY: Crater region X- and Y-coordinates.
%           inputs: Analysis input structure.
%           cutZ: Elevation grid cut to the edifice.
%           cutX: X-coordinate grid cut to the edifice.
%           cutY: Y-coordinate grid cut to the edifice.
%           Version: Analysis script version number.
%           StartTime: DateTime object giving the time the script was
%               started.
%           EndTime: DateTime object giving the time the script ended.
%       TopoParams:
%           ZHyps_Areas: Normalized area values for Z CDF.
%           ZHyps_Vals: Normalized Z CDF values.
%           Slope: Local slopes of DEM (GRIDobj).
%           Slope_Hyps_Areas: Normalized area values for slope CDF.
%           Slope_Hyps_Vals: Slope CDF values.
%           Curvature_Profile: Profile curvatures of DEM (GRIDobj).
%           CProf_Hyps_Areas: Normalized area values for profile curvature 
%               CDF.
%           CProf_Hyps_Vals: Profile curvature CDF values.
%           Curvature_Planform: Planform curvatures of DEM (GRIDobj).
%           CPlan_Hyps_Areas: Normalized area values for planform curvature 
%               CDF.
%           CPlan_Hyps_Vals: Profile planform CDF values.
%           Roughness: Indexed-structure of roughness analysis, includes:
%               Roughness_Grid = Grid of roughness values.
%               WindowSize = Pixel-size of roughness window.
%               ExpectedWindowRes = Expected, user-defined roughness
%                   window.
%               TrueWindowRes = Actual roughness based on grid resolution.
%               Hypsometry_Areas = Normalized area values for roughness.
%               Hypsometry_Values = Roughness CDF values.
%               Hypsometry_NormValues = Nomalized roughness CDF values.
%           TopographicCenter_XY: X-Y location of highest topography.
%           GeometricCenter_XY: X-Y location of center of boundary.
%           VolumetricCenter_XY: X-Y location of center of volume.
%       DrainageParams:
%           FD: Flow dirction object (FLOWobj).
%           A: Cumulative drainage area (GRIDobj).
%           D: Flow distance (GRIDobj).
%           DB: Separated drainage basins (GRIDobj).
%           DBxy: X-Y coordinates of the drainage basin boundaries, 
%               separated by NaNs.
%           DBxy_cell: X-Y coordinates of the drainage basin boundaries, 
%               separated as a cell array.
%           DB_Hyps_Topo: Normalized topography values for drainage basin CDF.
%           DB_Hyps_numDB: Drainage basin CDF values.
%           HackLawFit_BasinLength: Best-fit parameters for Hack's Law 
%               (relationship between drainage area and length) based on 
%               basin length.
%           HackLawFit_FlowLength: Best-fit parameters for Hack's Law 
%               (relationship between drainage area and length) based on 
%               basin flow length.
%           HackLawDeviation_BasinLength: Basin deviations from 
%               best-fitting Hack's Law parameters using basin length
%               (GRIDobj).
%           HackLawDeviation_FlowLength: Basin deviations from 
%               best-fitting Hack's Law parameters using basin flow length 
%               (GRIDobj).
%           Basin_TopN: Value used to define the top basins.
%           TopDBxy: X-Y coordinates of the largest drainage basin 
%               boundaries, separated by NaNs.
%           TopDBxy_cell: X-Y coordinates of the largest drainage basin 
%               boundaries, separated as a cell array.
%           TopN_basinIDs: Basin IDs of the top basins.
%           TopN_All_Area_Slope: Drainage area and slope values of the
%               largest N basins. Used for slope-area analysis.
%           TopN_Flow_Area_Slope: Drainage area and slope values of only
%               longest flow paths of the top N basins. Used for slope-area
%               analysis.                   
%           TopN_AreaThreshold_Rs: Array that contains the
%               dynamically-calculated drainage area threshold for
%               channels, the r^2 value of the hillslope regression, r^2
%               value of channel regression, and r' value used to determine
%               the threshold.
%           TopN_TransitionThreshold_A: Area threhold for transition zone
%               initiation.
%           TopN_MN: M/N value derived from linear regression of slope-area 
%               plots to determine best drainage area threshold for 
%               channelization.
%           Basin_Statistics: Matrix of basin IDs, drainage areas, max flow
%               lengths, basin lengths, widths, reliefs, orientations,
%               hypsometry integrals, mean slopes, and sinuosity.
%           Basin_Statistics_Grids: Structure of grids for overall basin
%               lengths, heights, widths, orientations, flow lengths,
%               hypsometry integrals, mean slopes, sinuosities, and 
%               drainage area.
%           Basin_Cross_Statistics: Array containing cross-basin stistics.
%               Rows are ordered as: basin ids, drainage sampling 
%               x-coordinates, drainage sampling y-coordinates, drainage
%               sampling elevations, cross-basin widths (perpendicular to
%               broad basin orientation) cross-basin relief (perpendicular
%               to broad basin orientation).
%           Basin_Statistics_Titles: Cell array of metrics titles
%               corresponding with Basin_Statistics columns.
%           Basin_Contour_ContourP_Count_Length_Area: Array containing 
%               analyzed contour values (in meters and percent relief from 
%               main flank), the number of basins in each contour, the
%               length of the contour, and the area of the contour.
%       ChannelParams:
%           ChannelThreshold: Drainage area threshold used for channel
%               distinction.
%           S: Determined river channels (STREAMobj).
%           DD: Channel drainage densities (needs to be paired with S).
%           MaxDD: Maximum basin drainage density (GRIDobj).
%           TotalDD: Cumulative drainage density of entire edifice.
%           Concavity_Streams: Streams used to analyze channel concavity.
%           Concavity_Stats: Stream concavity statistics, including
%               contributing drainage area, gradient, ks, and best-fit
%               concavity index (theta)
%           Concavity_BasinIDs: Basin IDs used for analysis.
%           Concavity_Mean: Mean concavity index of all analyzed streams.
%           Concavity_DEM: DEM of channel concavities (GRIDobj)
%           chiS: Determined river channels (STREAMobj) starting at the
%               same elevation, as required for Chi analysis.
%           chiS_topN: Determined river channels (STREAMobj) of the top N
%               basins, starting at the same elevation, as required for 
%               Chi analysis.
%           BestFit_MN: Best-fitting M/N ratio for Chi analysis. Either
%               provided by the user or derived from optimization.
%           Chi: Channel Chi values using the best-fitting m/n ratio
%               (needs to be paired with S).
%           MaxChi: Maximum basin Chi value (GRIDobj).
%           UpstreamChi: Chi values projected upstream to the divides.
%               (GRIDobj)
%       DivideParams:
%           DivideTopo: Divide catalog ordered by topography.
%               DVD: Divide structure (DIVIDEobj).
%               X: Divide X location.
%               Y: Divide Y location.
%           Divide_ReliefDistance: List of distances to common streams across
%               divides
%           Divide_AsymmetryIndex: Asymmetry index from Scherler & 
%               Schwanghart (2020)
%           Divide_ChiDifference: Differences in Chi across divides.
%           Junction_X_Y_C_Z_D_A: Divide junction information, including
%               X-coordinates, Y-coordinates, conjuction index, elevation,
%               distance along divide, and asymmetry index.
%           VerticalDistance: Grid of vertical distances between divides
%               and channels (GRIDobj).

%% Unpack input
pack = DrainageVolc_DefaultVals(pack);
disp(sprintf('USING DRAINAGEVOLC VERSION %s',pack.version))
disp('Unpacking Input and Importing Data...')

tifFile = pack.tifFile;
dx = pack.dx;
boundaryXY = pack.boundaryXY;
craterXY = pack.craterXY;
hypsIter = pack.hypsIter;
channelThreshold = pack.channelThreshold;
plotResults = pack.plotResults;
saveFigFolder = pack.saveFigFolder;
saveResFolder = pack.saveResFolder;
figPrefix = pack.figPrefix;
roughnessWindows = pack.roughnessWindows;
maskRegions = pack.maskMap;
basinTopN = pack.basinTopN;
parallelProc = pack.parallelProc;
dynThreshold = pack.dynamicThreshold;
dynThresholdStep = pack.dynamicThresholdPixelStep;
basinContIter = pack.basinContIter;
basinTopNFromContLength = pack.basinTopNFromContLength;
smoothBasinPointWavelength = pack.smoothBasinPointWavelength;
MN = pack.MN;

mets.GeneralParams.Version = pack.version;
mets.GeneralParams.StartTime = datetime('now');

%% Setup
% If boundary is given as shapefile, convert to an array.
if ischar(boundaryXY)
    Sh = shaperead(boundaryXY);
    tx = Sh.X;
    ty = Sh.Y;
    if isnan(tx(end))
        tx = tx(1:end-1);
        ty = ty(1:end-1);
    end
    
    try
        [ttx,tty,~] = ll2utm(ty,tx);
        boundaryXY = [ttx',tty'];
    catch
        warning('Unable to convert boundary from Lat/Lon, assuming already in UTM')
        boundaryXY = [tx',ty'];
    end
end

% If crater is given as shapefile, convert to cell array.
if ischar(craterXY)
    Sh = shaperead(craterXY);
    craterXY = [];
    cellCraterXY = {};
    for i = 1:size(Sh,1)
        tx = Sh(i).X;
        ty = Sh(i).Y;
        if isnan(tx(end))
            tx = tx(1:end-1);
            ty = ty(1:end-1);
        end

        try
            [ttx,tty,~] = ll2utm(ty,tx);
            craterXY = [craterXY;ttx',tty'];
            cellCraterXY = [cellCraterXY;{[ttx',tty']}];
        catch
            warning('Unable to convert crater from Lat/Lon, assuming already in UTM')
            craterXY = [craterXY;tx',ty'];
            cellCraterXY = [cellCraterXY;{[tx',ty']}];
        end
        craterXY = [craterXY;NaN,NaN];
    end
elseif ~isempty(craterXY)
    cellCraterXY = {};
    if iscell(craterXY)
        tmp = [];
        for i = 1:length(craterXY)
            tmp = [tmp;craterXY{i}];
            cellCraterXY = [cellCraterXY;{craterXY{i}}];
        end
        craterXY = tmp;
    else
        cellCraterXY = {craterXY};
    end
else
    cellCraterXY = {};
    craterXY = [];
end

% If mask is given as shapefile, convert to an array.
maskXY = [];
if ischar(maskRegions)
    Sh = shaperead(maskRegions);
    tx = Sh.X;
    ty = Sh.Y;
    if isnan(tx(end))
        tx = tx(1:end-1);
        ty = ty(1:end-1);
    end
    
    try
        [ttx,tty,~] = ll2utm(ty,tx);
        maskXY = [ttx',tty'];
    catch
        warning('Unable to convert mask from Lat/Lon, assuming already in UTM')
        maskXY = [tx',ty'];
    end
end

%% Load DEM
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

DEM0 = DEM;

% Fill sinks, except for given region
if isempty(craterXY)
    DEMf = fillsinks(DEM);
else
    [~,tx,ty] = GRIDobj2mat(DEM);
    [tX,tY] = meshgrid(tx,ty);
    tS = DEM;
    tS.Z = zeros(size(tS.Z));
    for i = length(cellCraterXY)
        p = inpolygon(tX,tY,cellCraterXY{i}(:,1),cellCraterXY{i}(:,2));
        tS.Z(p==1) = 1;
    end
    DEMf = fillsinks(DEM,tS);
end

[~,X,Y] = GRIDobj2mat(DEMf);
[XG,YG] = meshgrid(X,Y);
pp = inpolygon(XG,YG,boundaryXY(:,1),boundaryXY(:,2));
DEMf.Z(~pp) = NaN;

% If masking, change grid.
if ~isempty(maskXY)
    pp = inpolygon(XG,YG,maskXY(:,1),maskXY(:,2));
    DEMf.Z(pp) = NaN;
    
    pp = inpolygon(boundaryXY(:,1),boundaryXY(:,2),maskXY(:,1),maskXY(:,2));
    boundaryXY(pp,:) = [];
end

[Z,~,~] = GRIDobj2mat(DEMf);

mets.GeneralParams.DEM0 = DEM0;
mets.GeneralParams.X = X;
mets.GeneralParams.Y = Y;
mets.GeneralParams.Z = Z;
mets.GeneralParams.boundaryXY = boundaryXY;
mets.GeneralParams.craterXY = craterXY;
mets.GeneralParams.inputs = pack;

%% Separate Volcano
disp('Clipping Topography...')
% Cut DEM
DEMf = crop(DEMf);

[~,X,Y] = GRIDobj2mat(DEMf);
[XG,YG] = meshgrid(X,Y);
pp = inpolygon(XG,YG,boundaryXY(:,1),boundaryXY(:,2));
DEMf.Z(~pp) = NaN;
[Z,xc,yc] = GRIDobj2mat(DEMf);
[Xc,Yc] = meshgrid(xc,yc);
mets.GeneralParams.cutZ = Z;
mets.GeneralParams.cutX = X;
mets.GeneralParams.cutY = Y;

[ZHyps_Areas,ZHyps_Vals] = DrainageVolc_HypsometryValue(DEMf.Z,hypsIter,1);

mets.GeneralParams.DEM = DEMf;

mets.TopoParams.ZHyps_Areas = ZHyps_Areas;
mets.TopoParams.ZHyps_Vals = ZHyps_Vals;

%% Basic Topography Metrics
disp('Calculating Topography Metrics...')
% Slope
disp('   Slope...')
slp = gradient8(DEMf,'degree');
[Slope_Hyps_Areas,Slope_Hyps_Vals] = DrainageVolc_HypsometryValue(slp.Z,hypsIter,0);

mets.TopoParams.Slope = slp;
mets.TopoParams.Slope_Hyps_Areas = Slope_Hyps_Areas;
mets.TopoParams.Slope_Hyps_Vals = Slope_Hyps_Vals;   

% Curvature
disp('   Curvature...')
curve_prof = curvature(DEMf,'profc');
curve_plan = curvature(DEMf,'planc');

[CProf_Hyps_Areas,CProf_Hyps_Vals] = DrainageVolc_HypsometryValue(curve_prof.Z,hypsIter*.1,0);
[CPlan_Hyps_Areas,CPlan_Hyps_Vals] = DrainageVolc_HypsometryValue(curve_plan.Z,hypsIter*.1,0);

mets.TopoParams.Curvature_Profile = curve_prof;
mets.TopoParams.CProf_Hyps_Areas = CProf_Hyps_Areas;
mets.TopoParams.CProf_Hyps_Vals = CProf_Hyps_Vals;
mets.TopoParams.Curvature_Planform = curve_plan;
mets.TopoParams.CPlan_Hyps_Areas = CPlan_Hyps_Areas;
mets.TopoParams.CPlan_Hyps_Vals = CPlan_Hyps_Vals;

% Roughness
disp('   Roughness...')
roughnessVals = [];
for i = 1:length(roughnessWindows)
    roughI = round(roughnessWindows(i)/dx,0);
    if mod(roughI,2) == 0
        roughI = roughI+1;
    end
    roughGrid = roughness(DEMf,'tpi',[roughI,roughI]);
    [roughHypsAreas,roughHypsVals] = DrainageVolc_HypsometryValue(roughGrid.Z,hypsIter,0);
    [~,roughHypsNorm] = DrainageVolc_HypsometryValue(roughGrid.Z,hypsIter,1);
    
    roughnessVals(i).Roughness_Grid = roughGrid;
    roughnessVals(i).WindowSize = [roughI,roughI];
    roughnessVals(i).ExpectedWindowRes = roughnessWindows(i);
    roughnessVals(i).TrueWindowRes = roughI*dx;
    roughnessVals(i).Hypsometry_Areas = roughHypsAreas;
    roughnessVals(i).Hypsometry_Values = roughHypsVals;
    roughnessVals(i).Hypsometry_NormValues = roughHypsNorm;
end

mets.TopoParams.Roughness = roughnessVals;

% Topographic, Geometric, and Volumetric Centers

%   Topographic Center (peak)
disp('   Landform Centers...')
[ii,jj] = find(Z==max(Z(:)));
TC_X = XG(ii,jj);
TC_Y = YG(ii,jj);

%   Geometric Center (planform center of volcano)
%       The moments script calculates the volume-based mathmatical moments
%       of the elevation grid. Using 1 for all elevations within the
%       volcano removes the topographic weighting of the moments. This is
%       equivilant to taking the mean X and Y of all elevation coordinates.
tt = (Z>0)*1;
[~,~,~,~,~,com,~,~] = Calculate_Moments(X,Y',tt);
GC_X = com(1);
GC_Y = com(2);

%   Volumetric Center (volumetric center of mass, see Lerner et al. 2020 for
%       formal equation and description of use.
Z(isnan(Z)) = 0;
[~,~,~,~,~,com,~,~] = Calculate_Moments(X,Y',Z);
Z(Z==0) = NaN;
VC_X = com(1);
VC_Y = com(2);

mets.TopoParams.TopographicCenter_XY = [TC_X,TC_Y];
mets.TopoParams.GeometricCenter_XY = [GC_X,GC_Y];
mets.TopoParams.VolumetricCenter_XY = [VC_X,VC_Y];

%% Drainage Metrics
disp('Calculating Drainage Metrics...')
FD = FLOWobj(DEMf,'preprocess','none');

mets.DrainageParams.FD = FD;

% Drainage Area
disp('   Drainage Area...')
A = flowacc(FD);
A.Z(isnan(DEMf.Z)) = NaN;

% Flow Distance
disp('   Flow Distance...')
D = flowdistance(FD);
D.Z(isnan(DEMf.Z)) = NaN;

% Drainage Basin
disp('   Drainage Basins...')
DB = drainagebasins(FD);
DB = shufflelabel(DB);
DB.Z = single(DB.Z);
DB.Z(isnan(DEMf.Z)) = NaN;

[DBxy,DBxy_cell] = DrainageVolc_GetBasinBoundaries(DB,A,dx,channelThreshold);

[DB_Hyps_Topo,DB_Hyps_numDB] = DraingeVolc_Count_VolcTopoBasins(DEMf.Z,DB.Z,hypsIter);

mets.DrainageParams.A = A;
mets.DrainageParams.D = D;
mets.DrainageParams.DB = DB;
mets.DrainageParams.DBxy = DBxy;
mets.DrainageParams.DBxy_cell = DBxy_cell;
mets.DrainageParams.DB_Hyps_Topo = DB_Hyps_Topo;
mets.DrainageParams.DB_Hyps_numDB = DB_Hyps_numDB;

% Drainage Basin Stats
disp('   Drainage Basin Statistics...')
[tot_db_stats,tot_db_grids,cross_db_stats] = DrainageVolc_Collect_DrainageBasin_Stats_wCross(Z,slp,A,DB,FD,D,parallelProc,smoothBasinPointWavelength);

% Hack's Relationship (L = c*A^b) Using Basin Length
disp('   Hack''s Law Parameterization...')
tmpPF_Vals = [tot_db_stats(:,1),tot_db_stats(:,2),tot_db_stats(:,4)];
tmpPF_Vals(sum(isnan(tmpPF_Vals),2)>0,:) = [];
tmpPF_Vals(tmpPF_Vals(:,3)==0,:) = [];
powerFit = polyfit(log10(tmpPF_Vals(:,2)),log10(tmpPF_Vals(:,3)),1);

mets.DrainageParams.HackLawFit_BasinLength = [10^powerFit(2),powerFit(1)];

    % Length deviation from best-fitting power law
    lDev = log10(tmpPF_Vals(:,3)) - log10(10^powerFit(2)*tmpPF_Vals(:,2).^powerFit(1));
    lDev_Map = DB;
    lDev_Map.Z(:) = NaN;
    for i = 1:length(lDev)
        lDev_Map.Z(DB.Z==tmpPF_Vals(i,1)) = lDev(i);
    end
    
    mets.DrainageParams.HackLawDeviation_BasinLength = [tmpPF_Vals,lDev];
    mets.DrainageParams.HackLawDeviation_BasinLength_Map = lDev_Map;
    mets.DrainageParams.HackLawDeviation_Titles = {'Basin ID','Drainage Area','Length','Deviation'};

% Hack's Relationship (L = c*A^b) Using Flow Length
tmpPF_Vals = [tot_db_stats(:,1),tot_db_stats(:,2),tot_db_stats(:,3)];
tmpPF_Vals(sum(isnan(tmpPF_Vals),2)>0,:) = [];
tmpPF_Vals(tmpPF_Vals(:,3)==0,:) = [];
powerFit = polyfit(log10(tmpPF_Vals(:,2)),log10(tmpPF_Vals(:,3)),1);

mets.DrainageParams.HackLawFit_FlowLength = [10^powerFit(2),powerFit(1)];

    % Length deviation from best-fitting power law
    lDev = log10(tmpPF_Vals(:,3)) - log10(10^powerFit(2)*tmpPF_Vals(:,2).^powerFit(1));
    lDev_Map = DB;
    lDev_Map.Z(:) = NaN;
    for i = 1:length(lDev)
        lDev_Map.Z(DB.Z==tmpPF_Vals(i,1)) = lDev(i);
    end
    
    mets.DrainageParams.HackLawDeviation_FlowLength = [tmpPF_Vals,lDev];
    mets.DrainageParams.HackLawDeviation_FlowLength_Map = lDev_Map;

% Drainage Basins Per Contour
Basin_Contour_ContourP_Count_Length_Area = DrainageVolc_Determine_Basin_Per_Contour(DEMf,DB,basinContIter);
mets.DrainageParams.Basin_Contour_ContourP_Count_Length_Area = Basin_Contour_ContourP_Count_Length_Area;

% Drainage Area - Slope Plots
disp('   Slope-Area Relationship...')
if basinTopNFromContLength
    tmp = Basin_Contour_ContourP_Count_Length_Area;
    tmp(tmp(:,2) >= .9,:) = [];
    contLength = tmp(:,3)./tmp(:,4);
    diffcontLength = [0;diff(contLength)];

    ii = find(diffcontLength<0,1,'last');
    basinTopN = (tmp(ii,2))-1;
end
[topNAS,topNFlowAS,topDBxy,DBxy_cell,topN_basinIDs] = DrainageVolc_Collect_SlopeAreas(tot_db_stats,DB,D,FD,slp,A,dx,basinTopN,DEMf);
mets.DrainageParams.TopDBxy = topDBxy;
mets.DrainageParams.TopDBxy_cell = DBxy_cell;
mets.DrainageParams.TopN_All_Area_Slope = topNAS;
mets.DrainageParams.TopN_Flow_Area_Slope = topNFlowAS;
mets.DrainageParams.Basin_TopN = basinTopN;
mets.DrainageParams.TopN_basinIDs = sortrows(topN_basinIDs);

% Find Dynamic Threshold
if dynThreshold
    [basinAreaThreshold,transitionThreshold,MN] = DrainageVolc_SA_PiecewiseRegression(topNFlowAS,dynThresholdStep*dx^2);
    channelThreshold = basinAreaThreshold(1);
    
    mets.DrainageParams.TopN_AreaThreshold_Rs = basinAreaThreshold;
    mets.DrainageParams.TopN_TransitionThreshold_A = transitionThreshold;
    mets.DrainageParams.TopN_MN = abs(MN);
else
    mets.DrainageParams.TopN_AreaThreshold_Rs = [];
    mets.DrainageParams.TopN_TransitionThreshold_A = [];
    mets.DrainageParams.TopN_MN = [];
end
mets.ChannelParams.ChannelThreshold = channelThreshold;

mets.DrainageParams.Basin_Statistics = tot_db_stats;
mets.DrainageParams.Basin_Statistics_Grids = tot_db_grids;
mets.DrainageParams.Basin_Cross_Statistics = cross_db_stats;
mets.DrainageParams.Basin_Statistics_Titles = {'ID','Drainage Area','Flowpath Length','Basin Length','Basin Width','Basin Height','Basin Orientation','Basin Hypsometry Integral','Mean Basin Slope','Flowpath Sinuosity','Basin Euclidean Distance','Basin Length Distance to Largest Width'};

%% Channel Metrics
disp('Calculating Channel Metrics...')
S = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');

mets.ChannelParams.S = S;

% Drainage density
disp('   Drainage Density...')
DD = drainagedensity(S,FD,'nal');
DDG = drainagedensity(S,FD,'grid');
DB_DD = DB;
DB_DD.Z = zeros(size(DB.Z))*NaN;
DB_mDD = DB_DD;
uniDB = unique(DB.Z(:));
uniDB(isnan(uniDB)) = [];

sTmp = STREAMobj2GRIDobj(S);

for i = 1:length(uniDB)
    tmpDB = DB.Z == uniDB(i);
    maxVal = max(max(DDG.Z(tmpDB==1)));
    DB_mDD.Z(tmpDB==1) = maxVal;

    tmpZ = DEMf;
    tmpZ.Z(DB.Z~=uniDB(i)) = NaN;
    tmpZ.Z(~sTmp.Z) = NaN;
    [ii,jj] = find(tmpZ.Z==min(tmpZ.Z(:)),1);
    if ~isempty(ii)
        DB_DD.Z(tmpDB==1) = DDG.Z(ii,jj);
    end
end

mets.ChannelParams.DD = DD;
mets.ChannelParams.MaxDD = DB_mDD;
mets.ChannelParams.BasinDD = DB_DD;

mets.ChannelParams.TotalDD = DrainageVolc_Collect_Total_DrainageDensity(Z,S,A,DB,dx);

% Concavity Analysis
disp('   Channel Concavity...')
[topN_Ss,topN_S_Stats] = DrainageVolc_Collect_Stream_Concavity(tot_db_stats,DEMf,A,DB,FD,topN_basinIDs,channelThreshold);

mets.ChannelParams.Concavity_Streams = topN_Ss;
mets.ChannelParams.Concavity_Stats = topN_S_Stats;
mets.ChannelParams.Concavity_BasinIDs = topN_basinIDs;

conDEM = DEMf;
conDEM.Z(:,:) = NaN;
for i = 1:length(topN_basinIDs)
    conDEM.Z(DB.Z==topN_basinIDs(i)) = abs(topN_S_Stats(i).theta);
end

allTheta = [];
for i = 1:size(topN_Ss,1)
    allTheta = [allTheta;topN_S_Stats(i).theta];
end
mets.ChannelParams.Concavity_Mean = nanmean(allTheta);
mets.ChannelParams.Concavity_DEM = conDEM;

% Chi
%   This function uses Beysian optimization to find the best-fitting m/n
%   ratio for chi.
disp('   Chi Analysis...')

% Get streams at the same elevation (required for chi).
[tmpZ,tmpx,tmpy] = GRIDobj2mat(DEM0);
[tmpX,tmpY] = meshgrid(tmpx,tmpy);
boundaryZ = interp2(tmpX,tmpY,tmpZ,boundaryXY(:,1),boundaryXY(:,2));
chiDEM = DEMf;
chiDEM.Z(chiDEM.Z<max(boundaryZ)) = NaN;
for i = 1:length(cellCraterXY)
    p = inpolygon(XG,YG,cellCraterXY{i}(:,1),cellCraterXY{i}(:,2));
    chiDEM.Z(p==1) = NaN;
end
chiFD = FLOWobj(chiDEM,'preprocess','none');
chiS = STREAMobj(chiFD,'minarea',channelThreshold,'unit','mapunits');
mets.ChannelParams.chiS = chiS;

chiS_topN = klargestconncomps(chiS,length(topN_basinIDs));
mets.ChannelParams.chiS_topN = chiS_topN;

if ~isempty(chiS_topN.ix)
    try
        if isnan(MN)
            if ~plotResults
                try
                    [MN,~,~,~] = mnoptimvar(chiS_topN,DEMf,A,'varfun',@robustcov,'plot',false);
                catch
                    [MN,~,~,~] = mnoptimvar(chiS_topN,DEMf,A,'plot',false);
                end
            else
                figure('Name','M / N Optimization','units','normalized','outerposition',[0 0 1 1])
                try
                    [MN,~,~,~] = mnoptimvar(chiS_topN,DEMf,A,'varfun',@robustcov);
                catch
                    [MN,~,~,~] = mnoptimvar(chiS_topN,DEMf,A);
                end
                if ~isempty(saveFigFolder)
                    saveas(gcf,[saveFigFolder,figPrefix,'MN_Optimization.fig']);
                    saveas(gcf,[saveFigFolder,figPrefix,'MN_Optimization.png']);
                end
            end
        end
    
        Chi = chitransform(chiS,A,'a0',channelThreshold,'mn',MN);
        ChiG = STREAMobj2GRIDobj(chiS,Chi);
        DB_Chi = DB;
        DB_Chi.Z = zeros(size(DB.Z))*NaN;
        uniDB = unique(DB.Z(:));
        for i = 1:length(uniDB)
                tmpDB = DB.Z == uniDB(i);
                maxVal = max(max(ChiG.Z(tmpDB==1)));
                
                DB_Chi.Z(tmpDB==1) = maxVal;
        end

        Upstream_Chi = DrainageVolc_Project_Chi_Upstream(DEMf,DB,chiS,Chi);
    catch
        MN = NaN;
        Chi = NaN;
        DB_Chi = NaN;
        Upstream_Chi = NaN;
    end
else
    MN = NaN;
    Chi = NaN;
    DB_Chi = NaN;
    Upstream_Chi = NaN;
end

mets.ChannelParams.BestFit_MN = MN;
mets.ChannelParams.Chi = Chi;
mets.ChannelParams.MaxChi = DB_Chi;
mets.ChannelParams.UpstreamChi = Upstream_Chi;

%% Divide Metrics
disp('Calculating Divide Metrics...')
if ~pack.Analyze_Divides
    mets.DivideParams = [];
else

    % Collect divides
    DVD_Topo = DIVIDEobj(FD,S,'type','topo');
    DVD_Topo = cleanedges(DVD_Topo,FD);
    DVD_Topo = sort(DVD_Topo);
    
    % Order divides and find their roots
    [X_Topo,Y_Topo] = ind2coord(DVD_Topo,DVD_Topo.IX(1:end-1));
    
    Divide_Topo.DVD = DVD_Topo;
    Divide_Topo.X = X_Topo;
    Divide_Topo.Y = Y_Topo;
    mets.DivideParams.Divide_Topo = Divide_Topo;
    
    % Calculate Hillslope Asymmetry Values
    disp('   Hillslope Asymmetry...')
    vertDist = vertdistance2stream(FD,S,DEMf);
    vertDist.Z(isinf(vertDist.Z)) = NaN;
    
    [pp,qq] = getvalue(DVD_Topo,vertDist);
    reliefDistance = abs(diff(abs([pp,qq]),1,2));
    asymIndex = reliefDistance./sum([pp,qq],2);
    asymIndex(asymIndex<0) = 0;
    asymIndex(isinf(asymIndex)) = NaN;
    
    mets.DivideParams.Divide_ReliefDistance = reliefDistance;
    mets.DivideParams.Divide_AsymmetryIndex = asymIndex;

    % Calculate Chi differences across the divide
    disp('   Chi differences...')
    try
        [pp,qq] = getvalue(DVD_Topo,Upstream_Chi);
        chiDifference = abs(diff(abs([pp,qq]),1,2));
    catch
        chiDifference = NaN;
    end
    mets.DivideParams.Divide_ChiDifference = chiDifference;
    
    % Determine Divide Connectivity
    disp('   Divide Connectivity...')
    try
        [CJ,CJ_x,CJ_y] = jctcon(DVD_Topo);
        Junction_X_Y_C_Z_D_A = DrainageVolc_Collect_JunctionStats(DEMf,DVD_Topo,asymIndex,CJ,CJ_x,CJ_y);
    catch er
        Junction_X_Y_C_Z_D_A = [NaN,NaN,NaN,NaN,NaN,NaN];
    end
    
    mets.DivideParams.Junction_X_Y_C_Z_D_A = Junction_X_Y_C_Z_D_A;
    mets.DivideParams.VerticalDistance = vertDist;
end

%% Save All Results
mets.GeneralParams.EndTime = datetime('now');
if ~isempty(saveResFolder)
    disp('Saving Results...')
    save([saveResFolder,figPrefix,'DrainageVolc_Results.mat'],'mets')
end

%% Plot Results
if plotResults
    disp('Plotting Results...')
    DrainageVolc_Plots(mets);
end


end