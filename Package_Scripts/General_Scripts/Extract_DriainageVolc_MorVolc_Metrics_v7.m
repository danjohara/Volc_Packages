function [dMets,mMets] = Extract_DriainageVolc_MorVolc_Metrics_v7(drainageVolcFile,morVolcFile)
%%
% Name: Extract_DriainageVolc_MorVolc_Metrics
% Date: 08/20/2021 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to collect specific metrics from DrainageVolc and
%   MorVolc result files.
% 
% Input:
%   drainageVolcFile: Folder path and filename of DrainageVolc result file.
%   morVolcFile: Folder path and filename of MorVolc result file.
%
% Output:
%   dMets: Extracted DrainageVolc metrics.
%   mMets: Extracted MorVolc metrics.

%% Setup
warning('off','all')
% DrainageVolc Metrics
dMets.Ed_Area = NaN;
dMets.Ed_Area_km = NaN;
dMets.Basal_Width = NaN;
dMets.Max_Relief = NaN;

dMets.DB_Hyps_numDB = NaN;
dMets.DB_Hyps_Topo = NaN;
dMets.DB_Hyps_Integral = NaN;
dMets.Hacks_Law_Exponent_BasinLength = NaN;
dMets.Hacks_Law_Exponent_FlowLength = NaN;
dMets.Total_Drainage_Density = NaN;
dMets.Channel_Threshold = NaN;
dMets.Total_Conformity = NaN;

contStrings = {'00','10','20','30','40','50','60','70','80','90','95'};
contVals = [0:.1:.9,.95];
for i = 1:length(contVals)
    eval(sprintf('dMets.B_Count_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerEdArea_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerLength_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerArea_%s = NaN;',contStrings{i}));
end

radStrings = {'10','20','30','40','50','60','70','80','90','100'};
radVals = .1:.1:1;
for i = 1:length(radVals)
    eval(sprintf('dMets.RadB_Count_All_%s = NaN;',radStrings{i}));
    eval(sprintf('dMets.RadB_Count_Threshold_%s = NaN;',radStrings{i}));
end

dMets.Basin_IDs = NaN;
dMets.Basin_Lengths = NaN;
dMets.Basin_Reliefs = NaN;
dMets.Basin_Widths = NaN;
dMets.Basin_MaxWidthDists = NaN;
dMets.Basin_WidthAngles = NaN;
dMets.Basin_Slope = NaN;
dMets.Basin_Mean_Conformity = NaN;
dMets.Basin_Mean_ContSinuosity = NaN;
dMets.Basin_Hypsometry = NaN;
dMets.Basin_Slope_Variance = NaN;
dMets.Basin_Concavity = NaN;
dMets.Basin_MN = NaN;

dMets.Average_Basin_Length = NaN;
dMets.Average_Basin_Width = NaN;
dMets.Average_Basin_MaxWidthDists = NaN;
dMets.Average_Basin_WidthAngle = NaN;
dMets.Average_Basin_Relief = NaN;
dMets.Average_Basin_Slope = NaN;
dMets.Average_Basin_Conformity = NaN;
dMets.Average_Basin_Cont_Sinusity = NaN;
dMets.Average_Basin_Hypsometry = NaN;
dMets.Average_Basin_Slope_Variance = NaN;
dMets.Average_Basin_Concavity = NaN;
dMets.Average_Basin_MN = NaN;

dMets.P_Basin_IDs = NaN;
dMets.P_Basin_Lengths = NaN;
dMets.P_Basin_Reliefs = NaN;
dMets.P_Basin_Widths = NaN;
dMets.P_Basin_WidthAngles = NaN;
dMets.P_Basin_MaxWidthDists = NaN;
dMets.P_Basin_Slope = NaN;
dMets.P_Basin_Mean_Conformity = NaN;
dMets.P_Basin_Mean_ContSinuosity = NaN;
dMets.P_Basin_Hypsometry = NaN;
dMets.P_Basin_Slope_Variance = NaN;
dMets.P_Basin_Concavity = NaN;
dMets.P_Basin_MN = NaN;
dMets.P_Basin_Width_Distances = NaN;

dMets.P_Average_Basin_Length = NaN;
dMets.P_Average_Basin_Length_STD = NaN;
dMets.P_Average_Basin_Relief = NaN;
dMets.P_Average_Basin_Relief_STD = NaN;
dMets.P_Average_Basin_Width = NaN;
dMets.P_Average_Basin_Width_STD = NaN;
dMets.P_Average_Basin_Width_Angles = NaN;
dMets.P_Average_Basin_Width_Angles_STD = NaN;
dMets.P_Average_Basin_Hypsometry = NaN;
dMets.P_Average_Basin_Hypsometry_STD = NaN;
dMets.P_Average_Basin_Slope = NaN;
dMets.P_Average_Basin_Slope_STD = NaN;
dMets.P_Average_Basin_Mean_Conformity = NaN;
dMets.P_Average_Basin_Mean_Conformity_STD = NaN;
dMets.P_Average_Basin_Mean_ContSinuosity = NaN;
dMets.P_Average_Basin_Mean_ContSinuosity_STD = NaN;
dMets.P_Average_Basin_Concavity = NaN;
dMets.P_Average_Basin_Concavity_STD = NaN;
dMets.P_Average_Basin_MN = NaN;
dMets.P_Average_Basin_MN_STD = NaN;
dMets.P_Average_Basin_MN_CombinedStreams = NaN;

dMets.Norm_P_Basin_Lengths = NaN;
dMets.Norm_P_Basin_Reliefs = NaN;
dMets.Norm_P_Basin_Widths = NaN;

dMets.Norm_P_Average_Basin_Length = NaN;
dMets.Norm_P_Average_Basin_Relief = NaN;
dMets.Norm_P_Average_Basin_Width = NaN;

dMets.Norm_P_Average_Basin_Length_STD = NaN;
dMets.Norm_P_Average_Basin_Relief_STD = NaN;
dMets.Norm_P_Average_Basin_Width_STD = NaN;

dMets.P_All_Basin_Distances = NaN;
dMets.P_Average_Basin_Distance = NaN;
dMets.P_Average_Basin_Distance_STD = NaN;
dMets.P_Average_Basin_Distance_km = NaN;
dMets.P_Average_Basin_Distance_STD_km = NaN;
dMets.P_Basin_Spacing = NaN;

dMets.P_Average_Basin_Distance_Angle = NaN;
dMets.P_Average_Basin_Distance_Angle_STD = NaN;
dMets.P_Basin_Spacing_Angle = NaN;

dMets.Divide_Asymmetry_Gamma = NaN;

% MorVolc Metrics
mMets.Basal_Width = NaN;
mMets.Basal_Width_km = NaN;
mMets.Effective_Radius = NaN;
mMets.Effective_Radius_km = NaN;
mMets.Boundary_Perimeter = NaN;

mMets.Max_Radius = NaN;
mMets.Mean_Radius = NaN;
mMets.Max_Relief = NaN;
mMets.Max_Relief_km = NaN;

mMets.Hypsometry_Integral = NaN;
mMets.Total_Skewness = NaN;
mMets.Total_Kurtosis = NaN;

mMets.Topo_Vol_Center_Dist = NaN;
mMets.Topo_Geo_Center_Dist = NaN;
mMets.Vol_Geo_Center_Dist = NaN;

mMets.IDW_Volume = NaN;
mMets.Minimum_Eroded_Volume = NaN;

mMets.MainFlank_Irregularity_Indices_BFE = NaN;
mMets.MainFlank_Ellipticity_Indices_BFE = NaN;
mMets.MainFlank_Axis_Ellipticities = NaN;
mMets.MainFlank_Average_Irregularity_Indices_BFE = NaN;
mMets.MainFlank_Average_Ellipticity_Indices_BFE = NaN;
mMets.MainFlank_Average_Axis_Ellipticities = NaN;
mMets.MainFlank_Average_Irregularity_Indices_BFE_STD = NaN;
mMets.MainFlank_Average_Ellipticity_Indices_BFE_STD = NaN;
mMets.MainFlank_Average_Axis_Ellipticities_STD = NaN;

mMets.All_Average_Slope = NaN;
mMets.MainFlank_Average_Slope = NaN;
mMets.LowerFlank_Average_Slope = NaN;
mMets.Summit_Average_Slope = NaN;
mMets.All_STD_Slope = NaN;
mMets.MainFlank_STD_Slope = NaN;
mMets.LowerFlank_STD_Slope = NaN;
mMets.Summit_STD_Slope = NaN;

mMets.All_Slope_Variance = NaN;
mMets.MainFlank_Slope_Variance = NaN;
mMets.LowerFlank_Slope_Variance = NaN;
mMets.Summit_Slope_Variance = NaN;

MVcontStrings = {'10','20','30','40','50','60','70','80','90','95'};
MVcontVals = [.1:.1:.9,.95];
for i = 1:length(MVcontVals)
    eval(sprintf('mMets.Average_Slope_%s = NaN;',MVcontStrings{i}));
    eval(sprintf('mMets.Irregularity_BFE_%s = NaN;',MVcontStrings{i}));
    eval(sprintf('mMets.Ellipticity_BFE_%s = NaN;',MVcontStrings{i}));
    eval(sprintf('mMets.Axis_Ellipticity_%s = NaN;',MVcontStrings{i}));
end

mMets.Norm_Eroded_Volume = NaN;
mMets.Relief_Base_Ratio = NaN;


%% Load Data
try
    load(drainageVolcFile,'mets')
catch
    mets = [];
end

try
    load(morVolcFile,'morRes')
catch
    morRes = [];
end

%% Collect DrainageVolc Metrics
if isstruct(mets)
    %% Basic measurements
    dMets.Ed_Area = polyarea(mets.GeneralParams.boundaryXY(:,1),mets.GeneralParams.boundaryXY(:,2));
    dMets.Ed_Area_km = dMets.Ed_Area/1e6;
    dMets.Basal_Width = sqrt(dMets.Ed_Area/pi)*2;
    dMets.Max_Relief = range(mets.GeneralParams.cutZ(:));

    %% Basin topology metrics and distributions
    dMets.DB_Hyps_numDB = mets.DrainageParams.DB_Hyps_numDB;
    dMets.DB_Hyps_Topo = mets.DrainageParams.DB_Hyps_Topo;
    dMets.DB_Hyps_Integral = trapz(mets.DrainageParams.DB_Hyps_Topo,mets.DrainageParams.DB_Hyps_numDB./max(mets.DrainageParams.DB_Hyps_numDB));
    dMets.Hacks_Law_Exponent_BasinLength = mets.DrainageParams.HackLawFit_BasinLength(2);
    dMets.Hacks_Law_Exponent_FlowLength = mets.DrainageParams.HackLawFit_FlowLength(2);
    dMets.Total_Drainage_Density = mets.ChannelParams.TotalDD*1000;
    dMets.Channel_Threshold = mets.ChannelParams.ChannelThreshold./1e6;
    dMets.Total_Conformity = mets.ChannelParams.Conformity.Mean_Total_Conformity;

    %% Basin Topography Counts
    BCCpCLA = mets.DrainageParams.Basin_Contour_ContourP_Count_Length_Area;
    i0 = find(round(BCCpCLA(:,2),3)==0);
    for i = 1:length(contVals)
        tmp = BCCpCLA(:,2) - contVals(i);
        ii = find(round(abs(tmp),3) == 0,1);
        if ~isempty(ii)
            eval(sprintf('dMets.B_Count_%s = BCCpCLA(ii,3);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerEdArea_%s = BCCpCLA(ii,3)/(BCCpCLA(i0,5)/1e6);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerLength_%s = BCCpCLA(ii,3)/(BCCpCLA(ii,4)/1e3);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerArea_%s = BCCpCLA(ii,3)/(BCCpCLA(ii,5)/1e6);',contStrings{i}));
        end
    end

    %% Basin Radial Counts
    RA = mets.DrainageParams.Radial_Analysis.Basin_Count;
    for i = 1:length(radVals)
        tmp = RA(:,2) - radVals(i);
        ii = find(round(abs(tmp),3) == 0,1);
        if ~isempty(ii)
            eval(sprintf('dMets.RadB_Count_All_%s = RA(ii,3);',radStrings{i}));
            eval(sprintf('dMets.RadB_Count_Threshold_%s = RA(ii,4);',radStrings{i}));
        end
    end
    
    %% All Basin Data
    dMets.Basin_IDs = mets.DrainageParams.Basin_Statistics(:,1);
    dMets.Basin_Lengths = mets.DrainageParams.Basin_Statistics(:,4);
    dMets.Basin_Reliefs = mets.DrainageParams.Basin_Statistics(:,6);
    dMets.Basin_Widths = mets.DrainageParams.Basin_Statistics(:,5);
    dMets.Basin_MaxWidthDists = mets.DrainageParams.Basin_Statistics(:,12);
    dMets.Basin_WidthAngles = 2*atand(dMets.Basin_Widths./2./dMets.Basin_MaxWidthDists);
    dMets.Basin_Slope = mets.DrainageParams.Basin_Statistics(:,9);
    dMets.Basin_Hypsometry = mets.DrainageParams.Basin_Statistics(:,8);
    dMets.Basin_Slope_Variance = mets.DrainageParams.Basin_Statistics(:,13);

    dMets.Basin_Mean_Conformity = ones(size(dMets.Basin_IDs))*NaN;
    dMets.Basin_Mean_ContSinuosity = ones(size(dMets.Basin_IDs))*NaN;
    dMets.Basin_Concavity = ones(size(dMets.Basin_IDs))*NaN;
    dMets.Basin_MN = ones(size(dMets.Basin_IDs))*NaN;

    for i = 1:size(mets.ChannelParams.Conformity.Mean_Basin_Conformity,1)
        dMets.Basin_Mean_Conformity(dMets.Basin_IDs==mets.ChannelParams.Conformity.Mean_Basin_Conformity(i,1)) = mets.ChannelParams.Conformity.Mean_Basin_Conformity(i,2);
    end

    for i = 1:size(mets.DrainageParams.Basin_Contour_Sinuosity.BasinIDs,1)
        dMets.Basin_Mean_ContSinuosity(dMets.Basin_IDs==mets.DrainageParams.Basin_Contour_Sinuosity.BasinIDs(i)) = mets.DrainageParams.Basin_Contour_Sinuosity.Means(i);
    end

    for i = 1:size(mets.ChannelParams.Concavity_BasinIDs,1)
        dMets.Basin_Concavity(dMets.Basin_IDs==mets.ChannelParams.Concavity_BasinIDs(i)) = mets.ChannelParams.Concavity_Stats(i).theta;
    end

    for i = 1:length(mets.ChannelParams.Basin_Chi.Chi)
        dMets.Basin_MN(dMets.Basin_IDs==mets.ChannelParams.Basin_Chi.Chi(i).BasinID) = mets.ChannelParams.Basin_Chi.Chi(i).MN;
    end

    %% Mean of all basins
    dMets.Average_Basin_Length = nanmean(dMets.Basin_Lengths);
    dMets.Average_Basin_Width = nanmean(dMets.Basin_Widths);
    dMets.Average_Basin_MaxWidthDists = nanmean(dMets.Basin_MaxWidthDists);
    dMets.Average_Basin_WidthAngle = nanmean(dMets.Basin_WidthAngles);
    dMets.Average_Basin_Relief = nanmean(dMets.Basin_Reliefs);
    dMets.Average_Basin_Slope = nanmean(dMets.Basin_Slope);
    dMets.Average_Basin_Conformity = nanmean(dMets.Basin_Mean_Conformity);
    dMets.Average_Basin_Cont_Sinusity = nanmean(dMets.Basin_Mean_ContSinuosity);
    dMets.Average_Basin_Hypsometry = nanmean(dMets.Basin_Hypsometry);
    dMets.Average_Basin_Slope_Variance = nanmean(dMets.Basin_Slope_Variance);
    dMets.Average_Basin_Concavity = nanmean(dMets.Basin_Concavity);
    dMets.Average_Basin_MN = nanmean(dMets.Basin_MN);

    %% Top Basin Values
    dMets.P_Basin_IDs = mets.DrainageParams.TopN_basinIDs;
    tmp = zeros(size(dMets.P_Basin_IDs))*NaN;
    dMets.P_Basin_Lengths = tmp;
    dMets.P_Basin_Reliefs = tmp;
    dMets.P_Basin_Widths = tmp;
    dMets.P_Basin_WidthAngles = tmp;
    dMets.P_Basin_MaxWidthDists = tmp;
    dMets.P_Basin_Slope = tmp;
    dMets.P_Basin_Mean_Conformity = tmp;
    dMets.P_Basin_Mean_ContSinuosity = tmp;
    dMets.P_Basin_Hypsometry = tmp;
    dMets.P_Basin_Slope_Variance = tmp;
    dMets.P_Basin_Concavity = tmp;
    dMets.P_Basin_MN = tmp;
    dMets.P_Basin_Width_Distances = tmp;

    for i = 1:length(dMets.P_Basin_IDs)
        dMets.P_Basin_Lengths(i) = dMets.Basin_Lengths(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Reliefs(i) = dMets.Basin_Reliefs(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Widths(i) = dMets.Basin_Widths(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_WidthAngles(i) = dMets.Basin_WidthAngles(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_MaxWidthDists(i) = dMets.Basin_MaxWidthDists(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Slope(i) = dMets.Basin_Slope(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Mean_Conformity(i) = dMets.Basin_Mean_Conformity(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Mean_ContSinuosity(i) = dMets.Basin_Mean_ContSinuosity(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Hypsometry(i) = abs(dMets.Basin_Hypsometry(dMets.Basin_IDs == dMets.P_Basin_IDs(i)));
        dMets.P_Basin_Slope_Variance(i) = dMets.Basin_Slope_Variance(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_Concavity(i) = dMets.Basin_Concavity(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
        dMets.P_Basin_MN(i) = dMets.Basin_MN(dMets.Basin_IDs == dMets.P_Basin_IDs(i));
    end

    %% Average top basin values
    dMets.P_Average_Basin_Length = nanmean(dMets.P_Basin_Lengths);
    dMets.P_Average_Basin_Length_STD = nanstd(dMets.P_Basin_Lengths);
    dMets.P_Average_Basin_Relief = nanmean(dMets.P_Basin_Reliefs);
    dMets.P_Average_Basin_Relief_STD = nanstd(dMets.P_Basin_Reliefs);
    dMets.P_Average_Basin_Width = nanmean(dMets.P_Basin_Widths);
    dMets.P_Average_Basin_Width_STD = nanstd(dMets.P_Basin_Widths);
    dMets.P_Average_Basin_Width_Angles = nanmean(dMets.P_Basin_WidthAngles);
    dMets.P_Average_Basin_Width_Angles_STD = nanstd(dMets.P_Basin_WidthAngles);
    dMets.P_Average_Basin_Hypsometry = nanmean(dMets.P_Basin_Hypsometry);
    dMets.P_Average_Basin_Hypsometry_STD = nanstd(dMets.P_Basin_Hypsometry);
    dMets.P_Average_Basin_Slope = nanmean(dMets.P_Basin_Slope);
    dMets.P_Average_Basin_Slope_STD = nanstd(dMets.P_Basin_Slope);
    dMets.P_Average_Basin_Mean_Conformity = nanmean(dMets.P_Basin_Mean_Conformity);
    dMets.P_Average_Basin_Mean_Conformity_STD = nanstd(dMets.P_Basin_Mean_Conformity);
    dMets.P_Average_Basin_Mean_ContSinuosity = nanmean(dMets.P_Basin_Mean_ContSinuosity);
    dMets.P_Average_Basin_Mean_ContSinuosity_STD = nanstd(dMets.P_Basin_Mean_ContSinuosity);
    dMets.P_Average_Basin_Concavity = nanmean(dMets.P_Basin_Concavity);
    dMets.P_Average_Basin_Concavity_STD = nanstd(dMets.P_Basin_Concavity);
    dMets.P_Average_Basin_MN = nanmean(dMets.P_Basin_MN);
    dMets.P_Average_Basin_MN_STD = nanstd(dMets.P_Basin_MN);

    dMets.P_Average_Basin_MN_CombinedStreams = mets.ChannelParams.Total_Chi.MN;

    %% Normalized Basin Values
    dMets.Norm_P_Basin_Lengths = dMets.P_Basin_Lengths/(dMets.Basal_Width/2);
    dMets.Norm_P_Basin_Reliefs = dMets.P_Basin_Reliefs/dMets.Max_Relief;
    dMets.Norm_P_Basin_Widths = rad2deg(dMets.P_Basin_Widths/(dMets.Basal_Width/2*2*pi));

    dMets.Norm_P_Average_Basin_Length = nanmean(dMets.Norm_P_Basin_Lengths);
    dMets.Norm_P_Average_Basin_Relief = nanmean(dMets.Norm_P_Basin_Reliefs);
    dMets.Norm_P_Average_Basin_Width = nanmean(dMets.Norm_P_Basin_Widths);

    dMets.Norm_P_Average_Basin_Length_STD = nanstd(dMets.Norm_P_Basin_Lengths);
    dMets.Norm_P_Average_Basin_Relief_STD = nanstd(dMets.Norm_P_Basin_Reliefs);
    dMets.Norm_P_Average_Basin_Width_STD = nanstd(dMets.Norm_P_Basin_Widths);
    
    %% Basin spacing distances
    boundCumDists = 0;
    for i = 2:size(mets.GeneralParams.boundaryXY,1)
        tmpDist = sqrt((mets.GeneralParams.boundaryXY(i,1)-mets.GeneralParams.boundaryXY(i-1,1))^2 + (mets.GeneralParams.boundaryXY(i,2)-mets.GeneralParams.boundaryXY(i-1,2))^2);
        boundCumDists = [boundCumDists;boundCumDists(end)+tmpDist];
    end

    [DB,x,y] = GRIDobj2mat(mets.DrainageParams.DB);

    basinBoundDists = [];
    topBasinOutletPos = [];
    for i = 1:length(dMets.P_Basin_IDs)
        bb = bwboundaries(DB == dMets.P_Basin_IDs(i));
        tmpX = x(bb{1}(:,2));
        tmpY = y(bb{1}(:,1));
        tmpDists = pdist2([tmpX',tmpY],mets.GeneralParams.boundaryXY(:,1:2));
        [~,bI] = find(tmpDists == min(tmpDists(:)),1);
        topBasinOutletPos = [topBasinOutletPos;bI];
    end


    topBasinOutletPos = sortrows(topBasinOutletPos);
    for i = 2:length(topBasinOutletPos)
        basinBoundDists = [basinBoundDists;boundCumDists(topBasinOutletPos(i))-boundCumDists(topBasinOutletPos(i-1))];
    end
    basinBoundDists = [basinBoundDists;(boundCumDists(end)-boundCumDists(topBasinOutletPos(end))) + boundCumDists(topBasinOutletPos(1))];
    dMets.P_All_Basin_Distances = basinBoundDists;
    dMets.P_Average_Basin_Distance = nanmean(basinBoundDists);
    dMets.P_Average_Basin_Distance_STD = nanstd(basinBoundDists);
    
    dMets.P_Average_Basin_Distance_km = nanmean(basinBoundDists/1e3);
    dMets.P_Average_Basin_Distance_STD_km = nanstd(basinBoundDists/1e3);

    dMets.P_Basin_Spacing = (dMets.Basal_Width/2)/dMets.P_Average_Basin_Distance;

    %% Basin spacig angles
    dMets.P_Average_Basin_Distance_Angle = 360*mean(basinBoundDists)/(2*pi*dMets.Basal_Width/2);
    dMets.P_Average_Basin_Distance_Angle_STD = 360*std(basinBoundDists)/(2*pi*dMets.Basal_Width/2);
    dMets.P_Basin_Spacing_Angle = (dMets.Basal_Width/2)/dMets.P_Average_Basin_Distance_Angle;

    %% Divide Asymmetry PDF
    dMets.Divide_Asymmetry_Gamma = mets.DivideParams.DAI_Gamma.Gamma;
end

%% Collect MorVolc metrics
if isstruct(morRes)
    %% Basin shape 
    mMets.Basal_Width = sqrt(polyarea(morRes.GeneralParams.BoundaryXYZ(:,1),morRes.GeneralParams.BoundaryXYZ(:,2))/pi)*2;
    mMets.Basal_Width_km = mMets.Basal_Width/1e3;
    mMets.Effective_Radius = mMets.Basal_Width/2;
    mMets.Effective_Radius_km = mMets.Effective_Radius/1e3;
    mMets.Boundary_Perimeter = perimeter(polyshape(morRes.GeneralParams.BoundaryXYZ(:,1),morRes.GeneralParams.BoundaryXYZ(:,2)));

    rr = pdist2(morRes.TopoParams.TopographicCenter_XY,morRes.GeneralParams.BoundaryXYZ(:,1:2));
    mMets.Max_Radius = nanmax(rr);
    mMets.Mean_Radius = nanmean(rr);

    mMets.Max_Relief = range(morRes.GeneralParams.cutZ(:));
    mMets.Max_Relief_km = mMets.Max_Relief/1e3;

    %% Other shape parmeters
    mMets.Hypsometry_Integral = abs(trapz(morRes.TopoParams.ZHyps_Areas,morRes.TopoParams.ZHyps_Vals));
    mMets.Total_Skewness = sqrt(morRes.ShapeParams.Skewness(1)^2 + morRes.ShapeParams.Skewness(2)^2);
    mMets.Total_Kurtosis = sqrt(morRes.ShapeParams.Kurtosis(1)^2 + morRes.ShapeParams.Kurtosis(2)^2);

    mMets.Topo_Vol_Center_Dist = sqrt((morRes.TopoParams.TopographicCenter_XY(1)-morRes.TopoParams.VolumetricCenter_XY(1))^2 + ...
        (morRes.TopoParams.TopographicCenter_XY(2) - morRes.TopoParams.VolumetricCenter_XY(2))^2);
    mMets.Topo_Geo_Center_Dist = sqrt((morRes.TopoParams.TopographicCenter_XY(1)-morRes.TopoParams.GeometricCenter_XY(1))^2 + ...
        (morRes.TopoParams.TopographicCenter_XY(2) - morRes.TopoParams.GeometricCenter_XY(2))^2);
    mMets.Vol_Geo_Center_Dist = sqrt((morRes.TopoParams.VolumetricCenter_XY(1)-morRes.TopoParams.GeometricCenter_XY(1))^2 + ...
        (morRes.TopoParams.VolumetricCenter_XY(2) - morRes.TopoParams.GeometricCenter_XY(2))^2);

    %% Volumes
    if ~isempty(morRes.SizeParams.Volume.IDW)
        mMets.IDW_Volume = morRes.SizeParams.Volume.IDW;
    end
    
    mMets.Minimum_Eroded_Volume = morRes.SizeParams.Minimum_Eroded_Volume;

    %% Main Flank Contour shape metrics
    mMets.MainFlank_Irregularity_Indices_BFE = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_BFEllipse_Irregularity_Indexes];
    mMets.MainFlank_Ellipticity_Indices_BFE = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_BFEllipse_Ellipticity_Indexes];
    mMets.MainFlank_Axis_Ellipticities = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_Axis_Ellipticity];

    % Fix these indices to only look between first closed contour and
    % summit cutoff
    mMets.MainFlank_Irregularity_Indices_BFE(mMets.MainFlank_Irregularity_Indices_BFE(:,1)>morRes.GeneralParams.Summit_Contour,2) = NaN;
    mMets.MainFlank_Ellipticity_Indices_BFE(mMets.MainFlank_Ellipticity_Indices_BFE(:,1)>morRes.GeneralParams.Summit_Contour,2) = NaN;
    mMets.MainFlank_Axis_Ellipticities(mMets.MainFlank_Axis_Ellipticities(:,1)>morRes.GeneralParams.Summit_Contour,2) = NaN;
    
    mMets.MainFlank_Average_Irregularity_Indices_BFE = nanmean(mMets.MainFlank_Irregularity_Indices_BFE(:,2));
    mMets.MainFlank_Average_Ellipticity_Indices_BFE = nanmean(mMets.MainFlank_Ellipticity_Indices_BFE(:,2));
    mMets.MainFlank_Average_Axis_Ellipticities = nanmean(mMets.MainFlank_Axis_Ellipticities(:,2));

    mMets.MainFlank_Average_Irregularity_Indices_BFE_STD = nanstd(mMets.MainFlank_Irregularity_Indices_BFE(:,2));
    mMets.MainFlank_Average_Ellipticity_Indices_BFE_STD = nanstd(mMets.MainFlank_Ellipticity_Indices_BFE(:,2));
    mMets.MainFlank_Average_Axis_Ellipticities_STD = nanstd(mMets.MainFlank_Axis_Ellipticities(:,2));

    %% Slope Variance
    % All Edifice
    mMets.All_Average_Slope = morRes.SlopeParams.WholeEdifice_Mean_Slope;
    mMets.All_STD_Slope = morRes.SlopeParams.WholeEdifice_Std_Slope;
    mMets.All_Slope_Variance = mMets.All_STD_Slope/mMets.All_Average_Slope;

    % Main Flank
    mMets.MainFlank_Average_Slope = morRes.SlopeParams.Flank_Mean_Slope;
    mMets.MainFlank_STD_Slope = morRes.SlopeParams.Flank_Std_Slope;
    mMets.MainFlank_Slope_Variance = mMets.MainFlank_STD_Slope/mMets.MainFlank_Average_Slope;

    % Lower Flank
    mMets.LowerFlank_Average_Slope = morRes.SlopeParams.Lower_Flank_Mean_Slope;
    mMets.LowerFlank_STD_Slope = morRes.SlopeParams.Lower_Flank_Std_Slope;
    mMets.LowerFlank_Slope_Variance = mMets.LowerFlank_STD_Slope/mMets.LowerFlank_Average_Slope;

    % Summit
    mMets.Summit_Average_Slope = morRes.SlopeParams.Summit_Mean_Slope;
    mMets.Summit_STD_Slope = morRes.SlopeParams.Summit_Std_Slope;
    mMets.Summit_Slope_Variance = mMets.Summit_STD_Slope/mMets.Summit_Average_Slope;

    %% Contour Slopes, Irregularity, Ellipticity, and axis ellipticity
    [Z,~,~] = GRIDobj2mat(morRes.GeneralParams.DEM);
    Z = Z-min(Z(:));
    Z = Z./max(Z(:));

    for i = 2:length(MVcontVals)
        c1 = Convert_Contours(contourc(double(Z),[1,1]*MVcontVals(i)),1);
        if sqrt((c1(1,1)-c1(end,1))^2 + (c1(1,2)-c1(end,2))^2) < morRes.GeneralParams.DEM.cellsize
            ce = FitEllipse_Unwrap(c1);
            if isreal(ce.shortAxis)
                [ei_BFE,ii_BFE,~,~] = MorVolc_EI_II(c1,ce,1,1);
                ee = ce.shortAxis/ce.longAxis;
            else
                ei_BFE = NaN;
                ii_BFE = NaN;
                ee = NaN;
            end
        else
            ei_BFE = NaN;
            ii_BFE = NaN;
            ee = NaN;
        end

        SS = morRes.GeneralParams.S;
        SS(Z>=MVcontVals(i)) = NaN;
        SS(Z<MVcontVals(i-1)) = NaN;

        eval(sprintf('mMets.Average_Slope_%s = nanmean(SS(:));',MVcontStrings{i-1}));
        eval(sprintf('mMets.Irregularity_BFE_%s = ei_BFE;',MVcontStrings{i-1}));
        eval(sprintf('mMets.Ellipticity_BFE_%s = ii_BFE;',MVcontStrings{i-1}));
        eval(sprintf('mMets.Axis_Ellipticity_%s = ee;',MVcontStrings{i-1}));
    end

    %% Normalized Values
    mMets.Norm_Eroded_Volume = mMets.Minimum_Eroded_Volume/mMets.IDW_Volume;
    mMets.Relief_Base_Ratio = mMets.Max_Relief/(mMets.Effective_Radius);
end
warning('on','all')
end