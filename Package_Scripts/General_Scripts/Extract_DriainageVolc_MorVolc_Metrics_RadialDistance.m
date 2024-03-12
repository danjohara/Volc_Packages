function [dMets,mMets] = Extract_DriainageVolc_MorVolc_Metrics_RadialDistance(drainageVolcFile,morVolcFile,radDist,useMethod,reliefPer)
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
%   reliefPer: Relief percentile to define summit basins.
%
% Output:
%   dMets: Extracted DrainageVolc metrics.
%   mMets: Extracted MorVolc metrics.

%% Setup
% DrainageVolc Metrics
dMets.Hypsometry_Integral = NaN;
dMets.DB_Hyps_numDB = NaN;
dMets.DB_Hyps_Topo = NaN;
dMets.DB_Hyps_Integral = NaN;
dMets.Hacks_Law_Exponent_BasinLength = NaN;
dMets.Hacks_Law_Exponent_FlowLength = NaN;
dMets.Total_Drainage_Density = NaN;
dMets.Channel_Threshold = NaN;
dMets.Average_Basin_Length = NaN;
dMets.Average_Basin_Width = NaN;
dMets.Average_Basin_Relief = NaN;
dMets.Basin_Lengths = NaN;
dMets.Basin_Reliefs = NaN;
dMets.Basin_Widths = NaN;
dMets.Max_Radius = NaN;
dMets.Mean_Radius = NaN;

contStrings = {'00','10','20','30','40','50','60','70','80','90','95'};
contVals = [0:.1:.9,.95];
for i = 1:length(contVals)
    eval(sprintf('dMets.B_Count_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerEdArea_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerLength_%s = NaN;',contStrings{i}));
    eval(sprintf('dMets.B_CountPerArea_%s = NaN;',contStrings{i}));
end

radStrings = {'10','20','30','40','50','60','70','80','90','95'};
radVals = [.1:.1:.9,.95];
for i = 1:length(contVals)
    eval(sprintf('dMets.Rad_B_Count_%s = NaN;',contStrings{i}))
    eval(sprintf('dMets.Rad_B_CountPerLength_%s = NaN;',contStrings{i}))
end

dMets.P_Basin_Lengths = NaN;
dMets.P_Basin_Reliefs = NaN;
dMets.P_Basin_Widths = NaN;
dMets.P_Basin_Hypsometry = NaN;
dMets.P_Basin_Width_Distances = NaN;
dMets.P_Basin_Slope = NaN;

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
dMets.P_Average_Cross_Basin_Reliefs = NaN;
dMets.P_Average_Cross_Basin_Reliefs_STD = NaN;

dMets.P_All_Basin_Distances = NaN;
dMets.P_Average_Basin_Distance = NaN;
dMets.P_Average_Basin_Distance_STD = NaN;
dMets.P_Basin_Spacing = NaN;
dMets.P_Average_Basin_Distance_Angle = NaN;
dMets.P_Average_Basin_Distance_Angle = NaN;
dMets.P_Basin_Spacing_Angle = NaN;

dMets.P_Average_Basin_Distance_km = NaN;
dMets.P_Average_Basin_Distance_STD_km = NaN;
dMets.Ed_Area = NaN;
dMets.Ed_Area_km = NaN;

dMets.Divide_Asymmetry_Gamma = NaN;

% MorVolc Metrics
mMets.IDW_Volume = NaN;
mMets.Minimum_Eroded_Volume = NaN;
mMets.Max_Relief = NaN;
mMets.Boundary_Perimeter = NaN;
mMets.Basal_Width = NaN;
mMets.Irregularity_Indices_BFE = NaN;
mMets.Ellipticity_Indices_BFE = NaN;
mMets.Axis_Ellipticities = NaN;
mMets.Average_Irregularity_Indices_BFE = NaN;
mMets.Average_Ellipticity_Indices_BFE = NaN;
mMets.Average_Axis_Ellipticities = NaN;
mMets.Total_Skewness = NaN;

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

mMets.Effective_Radius = NaN;
mMets.Effective_Radius_km = NaN;

contStrings = {'10','20','30','40','50','60','70','80','90','95'};
contVals = [.1:.1:.9,.95];
for i = 1:length(contVals)
    eval(sprintf('mMets.Average_Slope_%s = NaN;',contStrings{i}));
    eval(sprintf('mMets.Irregularity_BFE_%s = NaN;',contStrings{i}));
    eval(sprintf('mMets.Ellipticity_BFE_%s = NaN;',contStrings{i}));
    eval(sprintf('mMets.Axis_Ellipticity_%s = NaN;',contStrings{i}));
end

% Combined, normalized metrics
dMets.Norm_P_Average_Basin_Length = NaN;
dMets.Norm_P_Average_Basin_Width = NaN;
dMets.Norm_P_Average_Basin_Relief = NaN;
dMets.Norm_P_STD_Basin_Relief = NaN;
dMets.Norm_P_STD_Basin_Length = NaN;
dMets.Norm_P_STD_Basin_Width = NaN;

mMets.Norm_Eroded_Volume = NaN;
mMets.Relief_Base_Ratio = NaN;

mMets.Basal_Width_km = NaN;
mMets.Max_Relief_km = NaN;


%% Load Data
try
    load(drainageVolcFile,'mets')
catch
    return;
end

try
    load(morVolcFile,'morRes')
catch
    return;
end

%% Collect DrainageVolc Metrics
if isstruct(mets)
    %% Topology metrics and distributions
    dMets.Hacks_Law_Exponent_BasinLength = mets.DrainageParams.HackLawFit_BasinLength(2);
    dMets.Hacks_Law_Exponent_FlowLength = mets.DrainageParams.HackLawFit_FlowLength(2);
    dMets.Total_Drainage_Density = mets.ChannelParams.TotalDD*1000;
    dMets.Basin_Lengths = mets.DrainageParams.Basin_Statistics(:,4);
    dMets.Basin_Reliefs = mets.DrainageParams.Basin_Statistics(:,6);
    dMets.Basin_Widths = mets.DrainageParams.Basin_Statistics(:,5);
    dMets.Channel_Threshold = mets.ChannelParams.ChannelThreshold./1e6;
    
    %% Mean topology values
    dMets.Average_Basin_Length = nanmean(dMets.Basin_Lengths);
    dMets.Average_Basin_Width = nanmean(dMets.Basin_Widths);
    dMets.Average_Basin_Relief = nanmean(dMets.Basin_Reliefs);
    
    %% Topography/basin distribution value
    dMets.Hypsometry_Integral = abs(trapz(mets.TopoParams.ZHyps_Areas,mets.TopoParams.ZHyps_Vals));
    dMets.DB_Hyps_numDB = mets.DrainageParams.DB_Hyps_numDB;
    dMets.DB_Hyps_Topo = mets.DrainageParams.DB_Hyps_Topo;
    dMets.DB_Hyps_Integral = trapz(mets.DrainageParams.DB_Hyps_Topo,mets.DrainageParams.DB_Hyps_numDB./max(mets.DrainageParams.DB_Hyps_numDB));
    dMets.Ed_Area = polyarea(mets.GeneralParams.boundaryXY(:,1),mets.GeneralParams.boundaryXY(:,2));
    
    %% Basin topology values and distributions for given percentile basins
    [Z,x,y] = GRIDobj2mat(mets.GeneralParams.DEM);
    [X,Y] = meshgrid(x,y);
    [distBins,BCounts_1,BCounts_2,~,nBCounts_1,nBCounts_2,BID_1,BID_2,~,~] = Calculate_BasinNum_ByDistance(mets,.01,'',0,0);
    ii = find(round(distBins,3) == round(radDist,3));
    if isempty(ii)
        error('radDist of %.2f not found in distance bins',radDist);
    end
    if useMethod==1
        DB_IDs = BID_1{ii};
        Rad_BCounts = BCounts_1;
        Rad_nBCounts = nBCounts_1;
    elseif useMethod == 2
        DB_IDs = BID_2{ii};
        Rad_BCounts = BCounts_2;
        Rad_nBCounts = nBCounts_2;
    end
    [DB,~,~] = GRIDobj2mat(mets.DrainageParams.DB);
    DB0 = DB;
    
%     Z = Z-min(Z(:));
%     Z = Z./max(Z(:));
%     DB(Z<reliefPer) = NaN;
%     DB_IDs = unique(DB(:));
%     DB_IDs(isnan(DB_IDs)) = [];
    
    dMets.P_Basin_Lengths = zeros(size(DB_IDs))*NaN;
    dMets.P_Basin_Reliefs = zeros(size(DB_IDs))*NaN;
    dMets.P_Basin_Widths = zeros(size(DB_IDs))*NaN;
    dMets.P_Basin_Hypsometry = zeros(size(DB_IDs))*NaN;
    dMets.P_Basin_Slope = zeros(size(DB_IDs))*NaN;
    dMets.P_Basin_Width_Distances = zeros(size(DB_IDs))*NaN;
    dMets.P_Cross_Basin_Reliefs = zeros(size(DB_IDs))*NaN;

    topBasinOutletPos = [];
    bxy = mets.GeneralParams.boundaryXY;
    newBXY = [];
    stepSize = .01;
    for j = 2:size(bxy,1)
        ddx = bxy(j,1)-bxy(j-1,1);
        ddy = bxy(j,2)-bxy(j-1,2);
        slp = ddy/ddx;
        hyp = sqrt(ddx^2+ddy^2);
        newHyps = (0:stepSize:1)*hyp;
    
        if ddx ~= 0
            newX = bxy(j-1,1)+sign(ddx)*newHyps*cos(atan(abs(slp)));
            newY = bxy(j-1,2)+sign(ddy)*newHyps*sin(atan(abs(slp)));
        else
            newY = bxy(j-1,2)+newHyps;
            newX = ones(size(newY))*bxy(j-1,1);
        end
        newBXY = [newBXY;newX',newY'];
    end
    bxy = newBXY;

    A = mets.DrainageParams.A;
    [Ag,~,~] = GRIDobj2mat(A);

    for i = 1:length(DB_IDs)
        ii = find(mets.DrainageParams.Basin_Statistics(:,1)==DB_IDs(i));
        dMets.P_Basin_Lengths(i) = mets.DrainageParams.Basin_Statistics(ii,4);
        dMets.P_Basin_Reliefs(i) = mets.DrainageParams.Basin_Statistics(ii,6);
        dMets.P_Basin_Widths(i) = mets.DrainageParams.Basin_Statistics(ii,5);
        dMets.P_Basin_Hypsometry(i) = mets.DrainageParams.Basin_Statistics(ii,8);
        dMets.P_Basin_Slope(i) = mets.DrainageParams.Basin_Statistics(ii,9);
        dMets.P_Basin_Width_Distances(i) = mets.DrainageParams.Basin_Statistics(ii,12);

        tmpA = Ag;
        tmpA(DB~=DB_IDs(i)) = NaN;
        [ii,jj] = find(tmpA==max(tmpA(:)),1);
        tmpDists = pdist2([X(ii,jj),Y(ii,jj)],bxy);
        [~,bI] = find(tmpDists == min(tmpDists(:)),1);
        topBasinOutletPos = [topBasinOutletPos;bI];

%         bb = bwboundaries(DB0 == DB_IDs(i));
%         tmpX = x(bb{1}(:,2));
%         tmpY = y(bb{1}(:,1));
%         tmpDists = pdist2([tmpX',tmpY],mets.GeneralParams.boundaryXY(:,1:2));
%         [~,bI] = find(tmpDists == min(tmpDists(:)),1);
%         topBasinOutletPos = [topBasinOutletPos;bI];
    end

    dMets.P_Average_Basin_Length = nanmean(dMets.P_Basin_Lengths);
    dMets.P_Average_Basin_Relief = nanmean(dMets.P_Basin_Reliefs);
    dMets.P_Average_Basin_Width = nanmean(dMets.P_Basin_Widths);
    dMets.P_Average_Basin_Width_Angles = nanmean(2*atand(dMets.P_Basin_Widths./2./dMets.P_Basin_Width_Distances)); % previously divided by P_Basin_Lengths (P_Basin_Width_Distances used in paper).
    dMets.P_Average_Basin_Hypsometry = abs(nanmean(dMets.P_Basin_Hypsometry));
    dMets.P_Average_Basin_Slope = nanmean(dMets.P_Basin_Slope);
    dMets.P_Average_Cross_Basin_Reliefs = nanmean(dMets.P_Cross_Basin_Reliefs);

    dMets.P_Average_Basin_Length_STD = nanstd(dMets.P_Basin_Slope);
    dMets.P_Average_Basin_Relief_STD = nanstd(dMets.P_Basin_Slope);
    dMets.P_Average_Basin_Width_STD = nanstd(dMets.P_Basin_Slope);
    dMets.P_Average_Basin_Width_Angles_STD = nanstd(2*atand(dMets.P_Basin_Widths./2./dMets.P_Basin_Width_Distances)); % previously divided by P_Basin_Lengths (P_Basin_Width_Distances used in paper).
    dMets.P_Average_Basin_Hypsometry_STD = nanstd(dMets.P_Basin_Hypsometry);
    dMets.P_Average_Basin_Slope_STD = nanstd(dMets.P_Basin_Slope);
    dMets.P_Average_Cross_Basin_Reliefs_STD = nanstd(dMets.P_Cross_Basin_Reliefs);
    
    %% Basin spacing distances
    boundCumDists = 0;
    for i = 2:size(bxy,1)
        tmpDist = sqrt((bxy(i,1)-bxy(i-1,1))^2 + (bxy(i,2)-bxy(i-1,2))^2);
        boundCumDists = [boundCumDists;boundCumDists(end)+tmpDist];
    end

    basinBoundDists = [];
    topBasinOutletPos = sortrows(topBasinOutletPos);
    for i = 2:length(topBasinOutletPos)
        basinBoundDists = [basinBoundDists;boundCumDists(topBasinOutletPos(i))-boundCumDists(topBasinOutletPos(i-1))];
    end
    basinBoundDists = [basinBoundDists;(boundCumDists(end)-boundCumDists(topBasinOutletPos(end))) + boundCumDists(topBasinOutletPos(1))];
    dMets.P_Average_Basin_Distance = mean(basinBoundDists);
    dMets.P_All_Basin_Distances = basinBoundDists;
    dMets.P_Average_Basin_Distance_STD = std(basinBoundDists);

    %% Basin count per contour length
    BCCpCLA = mets.DrainageParams.Basin_Contour_ContourP_Count_Length_Area;
    i0 = find(round(BCCpCLA(:,2),3)==0);
    
    contStrings = {'00','10','20','30','40','50','60','70','80','90','95'};
    contVals = [0:.1:.9,.95];
    for i = 1:length(contVals)
        tmp = BCCpCLA(:,2) - contVals(i);
        ii = find(round(abs(tmp),3) == 0,1);
%         ii = find(round(BCCpCLA(:,2),3)==contVals(i));
        if ~isempty(ii)
            eval(sprintf('dMets.B_Count_%s = BCCpCLA(ii,3);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerEdArea_%s = BCCpCLA(ii,3)/(BCCpCLA(i0,5)/1e6);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerLength_%s = BCCpCLA(ii,3)/(BCCpCLA(ii,4)/1e3);',contStrings{i}));
            eval(sprintf('dMets.B_CountPerArea_%s = BCCpCLA(ii,3)/(BCCpCLA(ii,5)/1e6);',contStrings{i}));
        end
    end

    %% Mean/max radii
    tc = mets.TopoParams.TopographicCenter_XY;
    if size(tc,2) > 2
        tc = [tc(1,1),tc(1,size(tc,2)/2+1)];
    end
    rr = pdist2(tc,mets.GeneralParams.boundaryXY(:,1:2));
    dMets.Max_Radius = nanmax(rr);
    dMets.Mean_Radius = nanmean(rr);

    %% Basin count per radial distance length
    radStrings = {'10','20','30','40','50','60','70','80','90','95'};
    radVals = [.1:.1:.9,.95];
    for i = 1:length(radVals)
        ii = find(round(distBins,3)==round(radVals(i),3));
        eval(sprintf('dMets.Rad_B_Count_%s = Rad_BCounts(ii);',radStrings{i}))
        eval(sprintf('dMets.Rad_B_CountPerLength_%s = Rad_nBCounts(ii);',radStrings{i}))
    end

    %% Divide Asymmetry PDF
    useOrder = 2;
    asBinWidths = .05;
    asBe = 0:asBinWidths:1;
    asBinMids = asBe(1:end-1)+asBinWidths/2;

    tmp = mets.DivideParams.Divide_AsymmetryIndex;
    to = divorder(mets.DivideParams.Divide_Topo.DVD);
    tmp(to.order < useOrder) = NaN;
    tmp(isnan(tmp)) = [];
    
    [n,~] = histcounts(tmp,asBe,'Normalization','pdf');
    dMets.Divide_Asymmetry_Gamma = trapz(asBinWidths,n);
end

%% Collect MorVolc metrics
if isstruct(morRes)
    %% Basic shape metrics
        % the summit region is erronous, so calculating everything based on
        % the height percentile
    [S,~,~] = GRIDobj2mat(mets.TopoParams.Slope);
    [Z,~,~] = GRIDobj2mat(mets.GeneralParams.DEM);

%     SS = S;
%     SS(Z>heightCutoff) = NaN;
%     SS(Z<max(morRes.GeneralParams.BoundaryXYZ(:,3))) = NaN;
%     mMets.Main_Flank_Average_Slope = nanmean(SS(:));
%     mMets.Main_Flank_STD_Slope = nanstd(SS(:));

    heightCutoff = range(Z(:))*reliefPer+min(Z(:));
    if ~isempty(morRes.SizeParams.Volume.IDW)
        mMets.IDW_Volume = morRes.SizeParams.Volume.IDW;
    end
    
    mMets.Minimum_Eroded_Volume = morRes.SizeParams.Minimum_Eroded_Volume;
    mMets.Max_Relief = morRes.SizeParams.Maximum_Height;
    mMets.Basal_Width = morRes.SizeParams.Basal_Width;
    mMets.Total_Skewness = sqrt(morRes.ShapeParams.Skewness(1)^2 + morRes.ShapeParams.Skewness(2)^2);
    
    pgon=polyshape(morRes.GeneralParams.BoundaryXYZ(:,1),morRes.GeneralParams.BoundaryXYZ(:,2));
    mMets.Boundary_Perimeter = perimeter(pgon);

    %% Contour shape metrics
    mMets.Irregularity_Indices_BFE = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_BFEllipse_Irregularity_Indexes];
    mMets.Ellipticity_Indices_BFE = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_BFEllipse_Ellipticity_Indexes];
    mMets.Axis_Ellipticities = [morRes.ShapeParams.Contour_Values,morRes.ShapeParams.Contour_Axis_Ellipticity];

    % Fix these indices to only look between first closed contour and
    % summit cutoff
    mMets.Irregularity_Indices_BFE(mMets.Irregularity_Indices_BFE(:,1)>heightCutoff,2) = NaN;
    mMets.Ellipticity_Indices_BFE(mMets.Ellipticity_Indices_BFE(:,1)>heightCutoff,2) = NaN;
    mMets.Axis_Ellipticities(mMets.Axis_Ellipticities(:,1)>heightCutoff,2) = NaN;
    
    mMets.Average_Irregularity_Indices_BFE = nanmean(mMets.Irregularity_Indices_BFE(:,2));
    mMets.Average_Ellipticity_Indices_BFE = nanmean(mMets.Ellipticity_Indices_BFE(:,2));
    mMets.Average_Axis_Ellipticities = nanmean(mMets.Axis_Ellipticities(:,2));

    %% Slope Variance
    % All Edifice
    mMets.All_Average_Slope = nanmean(S(:));
    mMets.All_STD_Slope = nanstd(S(:));
    mMets.All_Slope_Variance = mMets.All_STD_Slope/mMets.All_Average_Slope;

    % Main Flank
    SS = S;
    SS(Z>heightCutoff) = NaN;
    SS(Z<max(morRes.GeneralParams.BoundaryXYZ(:,3))) = NaN;
    mMets.MainFlank_Average_Slope = nanmean(SS(:));
    mMets.MainFlank_STD_Slope = nanstd(SS(:));
    mMets.MainFlank_Slope_Variance = mMets.MainFlank_STD_Slope/mMets.MainFlank_Average_Slope;

    % Lower Flank
    SS = S;
    SS(Z>max(morRes.GeneralParams.BoundaryXYZ(:,3))) = NaN;
    mMets.LowerFlank_Average_Slope = nanmean(SS(:));
    mMets.LowerFlank_STD_Slope = nanstd(SS(:));
    mMets.LowerFlank_Slope_Variance = mMets.LowerFlank_STD_Slope/mMets.LowerFlank_Average_Slope;

    % Summit
    SS = S;
    SS(Z<heightCutoff) = NaN;
    mMets.Summit_Average_Slope = nanmean(SS(:));
    mMets.Summit_STD_Slope = nanstd(SS(:));
    mMets.Summit_Slope_Variance = mMets.Summit_STD_Slope/mMets.Summit_Average_Slope;

%     mMets.All_Slope_Variance = morRes.SlopeParams.WholeEdifice_Std_Slope/morRes.SlopeParams.WholeEdifice_Mean_Slope;
%     try
%         mMets.MainFlank_Slope_Variance = morRes.SlopeParams.Flank_Std_Slope/morRes.SlopeParams.Flank_Mean_Slope;
%     catch
%         mMets.MainFlank_Slope_Variance = morRes.SlopeParams.Flank_Slp_Slope/morRes.SlopeParams.Flank_Mean_Slope;
%     end
%     mMets.LowerFlank_Slope_Variance = morRes.SlopeParams.Lower_Flank_Std_Slope/morRes.SlopeParams.Lower_Flank_Mean_Slope;
%     mMets.Summit_Slope_Variance = morRes.SlopeParams.Summit_Std_Slope/morRes.SlopeParams.Summit_Mean_Slope;

    %% Contour Slopes, Irregularity, Ellipticity, and axis ellipticity
    [Z,~,~] = GRIDobj2mat(mets.GeneralParams.DEM);
    Z = Z-min(Z(:));
    Z = Z./max(Z(:));

    contStrings = {'10','20','30','40','50','60','70','80','90','95'};
    contVals = [0:.1:.9,.95];
    for i = 2:length(contVals)
        c1 = Convert_Contours(contourc(double(Z),[1,1]*contVals(i)),1);
        if sqrt((c1(1,1)-c1(end,1))^2 + (c1(1,2)-c1(end,2))^2) < mets.GeneralParams.DEM.cellsize
            ce = FitEllipse_Unwrap(c1);
            if isreal(ce.shortAxis)
                [ei_BFE,ii_BFE,~,~] = Morvolc_Ellipticity_Irregularity_BFEllipse(c1(:,1),c1(:,2),ce);
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

        SS = S;
        SS(Z>=contVals(i)) = NaN;
        SS(Z<contVals(i-1)) = NaN;

        eval(sprintf('mMets.Average_Slope_%s = nanmean(SS(:));',contStrings{i-1}));
        eval(sprintf('mMets.Irregularity_BFE_%s = ei_BFE;',contStrings{i-1}));
        eval(sprintf('mMets.Ellipticity_BFE_%s = ii_BFE;',contStrings{i-1}));
        eval(sprintf('mMets.Axis_Ellipticity_%s = ee;',contStrings{i-1}));
    end
end

%% Collect normalized metrics
dMets.Norm_P_Basin_Lengths = dMets.P_Basin_Lengths/(mMets.Basal_Width/2);
dMets.Norm_P_Basin_Reliefs = dMets.P_Basin_Reliefs/mMets.Max_Relief;
dMets.Norm_P_Basin_Widths = rad2deg(dMets.P_Basin_Widths/(sqrt(dMets.Ed_Area/pi)*2*pi));
% dMets.Norm_P_Basin_Widths = dMets.P_Basin_Widths/P_Basin_Lengths;
% % dMets.Norm_P_Basin_Widths = 360*(dMets.P_Basin_Widths/mMets.Boundary_Perimeter);
% % dMets.Norm_P_Basin_Widths = rad2deg(dMets.P_Basin_Widths/(dMets.Mean_Radius*2*pi));

dMets.Norm_P_Average_Basin_Relief = nanmean(dMets.Norm_P_Basin_Reliefs);
dMets.Norm_P_Average_Basin_Length = nanmean(dMets.Norm_P_Basin_Lengths);
dMets.Norm_P_Average_Basin_Width = nanmean(dMets.Norm_P_Basin_Widths);

dMets.Norm_P_STD_Basin_Relief = nanstd(dMets.Norm_P_Basin_Reliefs);
dMets.Norm_P_STD_Basin_Length = nanstd(dMets.Norm_P_Basin_Lengths);
dMets.Norm_P_STD_Basin_Width = nanstd(dMets.Norm_P_Basin_Widths);

dMets.P_Basin_Spacing = (mMets.Basal_Width/2)/dMets.P_Average_Basin_Distance;

mMets.Norm_Eroded_Volume = mMets.Minimum_Eroded_Volume/(mMets.IDW_Volume+mMets.Minimum_Eroded_Volume)*100;
mMets.Relief_Base_Ratio = mMets.Max_Relief/(mMets.Basal_Width/2);

dMets.P_Average_Basin_Distance_Angle = 360*mean(basinBoundDists)/(2*pi*mMets.Basal_Width/2);
dMets.P_Average_Basin_Distance_Angle_STD = 360*std(basinBoundDists)/(2*pi*mMets.Basal_Width/2);
dMets.P_Basin_Spacing_Angle = (mMets.Basal_Width/2)/dMets.P_Average_Basin_Distance_Angle;

mMets.Basal_Width_km = mMets.Basal_Width/1e3;
mMets.Max_Relief_km = mMets.Max_Relief/1e3;
mMets.Effective_Radius = mMets.Basal_Width/2;
mMets.Effective_Radius_km = mMets.Effective_Radius/1e3;

dMets.P_Average_Basin_Distance_km = dMets.P_Average_Basin_Distance/1e3;
dMets.P_Average_Basin_Distance_STD_km = dMets.P_Average_Basin_Distance_STD/1e3;
dMets.Ed_Area_km = dMets.Ed_Area/1e6;
end