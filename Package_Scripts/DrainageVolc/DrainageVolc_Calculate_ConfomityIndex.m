function [totalMeanConformity,basinMeanConformity,DB_Conformity_Map,basinOrderMeanConformity,All_Segment_Info,DEMf] = DrainageVolc_Calculate_ConfomityIndex(DEM0,DEM,boundaryXY,DB,channelThreshold,wLength,sDist,verbose)
% Name: DrainageVolc_Calculate_ConfomityIndex
% Author: Daniel O'Hara
% Date: 06/09/2024 (mm/dd/yyyy)
% Description: Script to calculate conformity indices (Black et al., 2017)
%   over varying scales. This includes the mean conformity index of the
%   entire edifice, mean conformity indices of every applicable basin, and
%   mean conformity index of every stream order within a basin. The
%   conformity index measures the difference in orientation between river
%   channels and the path flow takes over a filtered topography, given as 
%   the positive cosine of the difference betweeen orientaitons. River 
%   channel and flow path orientations are calculated over the distance
%   given by sDist. Topography is filtered over the wavelength given by 
%   wLength; if it is NaN, wLength is set as the effective radius (radius
%   of a circle that has the same area as the edifice boundary).
%
%   This script uses parallel looping.
%
% Input:
%   DEM0: Initial GRIDobj of DEM. This is the topography that is filtered.
%   DEM: Cut GRIDobj of DEM.
%   boundaryXY: Array of boundary x-y coordinates.
%   DB: GRIDobj of drainage basins.
%   channelThreshold: Drainage area threshold for channelization.
%   wLength: Wavelength to filter topography. If NaN, wLength is calculated 
%       as the effective edifice radius.   
%   sDist: Channel and flow path distance over which orientation is
%       calculated.
%   verbose: Flag for whether updates should be outputted.
%
% Output:
%   totalMeanConformity: Mean conformity index of the entire edifice.
%   basinMeanConformity: Mean conformity indices of each basin with a
%       channel. First column is basin ID, second column is conformity.
%   DB_Conformity_Map: GRIDobj of basins, filled in by mean conformity
%       index.
%   basinOrderMeanConformity: Structure conformity indices for each basin
%       and Strahler order. Contains basin IDs, arrays of stream order
%       numbers, and matching arrays of conformity index.
%   All_Segment_Info: Structure of all channel segments used to determine
%       conformity indices. Contains basin IDs, segment XYZ and distance 
%       values, segment Strahler order, segment orientation, associated 
%       long-wavelength flowpath orientation, orientation differences, and
%       conformity index.
%   DEMf: Filtered topography.

warning('off','all')
%% Determine filter wavelength
if isnan(wLength)
    pa = polyarea(boundaryXY(:,1),boundaryXY(:,2));
    wLength = sqrt(pa/pi);
end

%% Filter topography
kern = ceil(wLength/DEM0.cellsize);
if mod(kern,2) == 0
    kern = kern+1;
end
DEM0f = filter(DEM0,'mean',[1,1]*kern);
[Z0f,X0,Y0] = GRIDobj2mat(DEM0f);
[X0,Y0] = meshgrid(X0,Y0);
si = scatteredInterpolant(double(X0(:)),double(Y0(:)),double(Z0f(:)),'natural');

[~,X,Y] = GRIDobj2mat(DEM);
[X,Y] = meshgrid(X,Y);
Zf = si(X,Y);
DEMf = GRIDobj(X,Y,Zf);
DEMf.Z(isnan(DEM.Z)) = NaN;

FDf = FLOWobj(DEMf);

%% Run through basins collecting conformity indices
uniDB = unique(DB.Z(:));
uniDB(isnan(uniDB)) = [];

DB_Conformity_Map = DB;
DB_Conformity_Map.Z(:) = NaN;

basinMeanConformity = [];
basinOrderMeanConformity = [];

All_Segment_Info = []; % will contain basin ID, x-y values, order, orientation, filtered orientation, orientation difference (gamma), conformity (cos(gamma))
allSegcounter = 1;
basinCounter = 1;

allConformities = [];

curPer = 0;
for i = 1:length(uniDB)
    if verbose > 0 && i/length(uniDB)>=curPer
        disp(sprintf('         %d%% Complete (%d / %d)',round(curPer*100,0),i,length(uniDB)))
        curPer = curPer + .1;
    end
    
    %% Test and Isolate basin
    tmpDEM = DEM;
    tmpDEM.Z(DB.Z~=uniDB(i)) = NaN;

    if sum(~isnan(tmpDEM.Z(:)))*tmpDEM.cellsize^2 < channelThreshold
        continue;
    end

    tmpDEM = crop(tmpDEM);

    %% Get basin STREAMobj
    FD = FLOWobj(tmpDEM);
    S = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');
    SO = streamorder(S);

    if isempty(S.ix)
        continue;
    end

    all_basinScale_confomities = [];
    all_basinScale_SO_conformities = cell(max(SO),1);
    for j = 1:max(SO)
        all_basinScale_SO_conformities{j} = [];
    end

    %% Seperate by tributary
    Ss = split(S);
    CS_m = STREAMobj2cell(Ss);

    for j = 1:length(CS_m)
        %% Separate by segment
        CS_2 = STREAMobj2cell(CS_m{j},'segment',sDist);

        for k = 1:length(CS_2)
            %% Collect segment order
            CS = CS_2{k};
            t1 = S.x == CS.x(1);
            t2 = S.y == CS.y(1);
            tt = t1.*t2;
            ii = find(tt == 1,1);
            CS_SO = SO(ii);

            %% Collect segment orientation
            tmpZ = getnal(CS,tmpDEM);
            CS_XYZ_dist = sortrows([CS.x,CS.y,tmpZ,CS.distance],3);

            CS_Phi = Calculate_Azimuths(CS_XYZ_dist(end,1:2),CS_XYZ_dist(1,1:2)); % CS is ordered from base level.

            %% Determine filtered topography flowpath and orientation
            R = sqrt((X-CS_XYZ_dist(end,1)).^2 + (Y-CS_XYZ_dist(end,2)).^2);
            ii = find(R == min(R(:)),1);
            [~,fp_dists,fp_x,fp_y] = flowpathextract(FDf,ii);
            fp_xy_dists = [fp_x,fp_y,fp_dists];
            fp_xy_dists(fp_xy_dists(:,3)>CS_XYZ_dist(end,4),:) = [];

            fp_endPoints = [fp_xy_dists(end,:);fp_xy_dists(1,:)];

            fp_Phi = Calculate_Azimuths(fp_xy_dists(1,1:2),fp_xy_dists(end,1:2)); % fp is ordered from flow path head, so need to reverse order of subtraction.
            
            %% Combine Data into array
            conformity = abs(cosd(CS_Phi-fp_Phi));
            All_Segment_Info(allSegcounter).Basin_ID = uniDB(i);
            All_Segment_Info(allSegcounter).X = CS_XYZ_dist(:,1);
            All_Segment_Info(allSegcounter).Y = CS_XYZ_dist(:,2);
            All_Segment_Info(allSegcounter).Z = CS_XYZ_dist(:,3);
            All_Segment_Info(allSegcounter).Dist = CS_XYZ_dist(:,4);
            All_Segment_Info(allSegcounter).Order = CS_SO;
            All_Segment_Info(allSegcounter).Orientation = CS_Phi;
            All_Segment_Info(allSegcounter).LongWavelength_EndPoints = fp_endPoints;
            All_Segment_Info(allSegcounter).LongWavelength_Orientation = fp_Phi;
            All_Segment_Info(allSegcounter).Orientation_Difference = abs(CS_Phi-fp_Phi);
            All_Segment_Info(allSegcounter).Conformity = conformity;

            allConformities = [allConformities;conformity];
            all_basinScale_confomities = [all_basinScale_confomities;conformity];
            all_basinScale_SO_conformities{CS_SO} = [all_basinScale_SO_conformities{CS_SO};conformity];

            allSegcounter = allSegcounter + 1;
        end
    end
    basinMeanConformity = [basinMeanConformity;uniDB(i),nanmean(all_basinScale_confomities)];
    bOrders = [1:max(SO)]';
    meanBOCon = zeros(max(SO),1);
    for j = 1:max(SO)
        manBOCon(j) = nanmean(all_basinScale_SO_conformities{j});
    end
    basinOrderMeanConformity(basinCounter).Basin_ID = uniDB(i);
    basinOrderMeanConformity(basinCounter).Stream_Orders = bOrders;
    basinOrderMeanConformity(basinCounter).Mean_Conformities = manBOCon;

    DB_Conformity_Map.Z(DB.Z==uniDB(i)) = nanmean(all_basinScale_confomities);

    basinCounter = basinCounter+1;
end

totalMeanConformity = nanmean(allConformities);
warning('on','all')
end