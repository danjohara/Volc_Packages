function [tot_db_stats,tot_db_grids,cross_db_stats] = DrainageVolc_Collect_DrainageBasin_Stats_wCross(Z,Slp,A,DB,FD,D,usePar,smoothBasinPointWavelength)
% Name: Collect_DrainageBasin_Stats_wCross
% Author: Daniel O'Hara
% Date: 03/03/2021 (mm/dd/yyyy)
% Description: Script to determine geometry statistics for drainage basins.
%   Geometries are calculated by considering the trendline between the
%   lowest point (P1) of the basin and the furthest location away (P2; 
%   generally corresponding with the highest point). Height is calculated 
%   as the range of elevations within the basin; length is midpoint line 
%   distance between P1 and P2; phi is the azimuthal angle from P2 to P1 
%   (e.g., down-slope map angle); and width is calculated by first rotating 
%   the data to set phi horizontal (e.g., 0 degrees in unit geometry), and 
%   calculting the range of bins of points between P1 and P2. Cross-statics 
%   (height and width) are calculated along the longest flow network. 
%   Values represent the basin height and width perpendicular to the 
%   general drainage basin orientation.
%
%   This script uses parallel looping.
%
% Input:
%   Z: Grid of elevations.
%   Slp: GRIDobj of slopes.
%   A: GRIDobj of upstream drainage areas.
%   DB: GRIDobj of drainage basins.
%   FD: FLOWobj of flow directions.
%   D: GRIDobj of flow distances.
%   usePar: Flag for whether parallel processing should be used.    
%   smoothBasinPointWavelength: Wavelength to smooth mid-basin points,
%       to help correct noisy data that can skew basin length values.
%
% Output:
%   tot_db_stats: Matrix of basin geometry statistics, basin IDs, drainage 
%       areas, max flow lengths, basin lengths, widths, reliefs, 
%       orientations, hypsometry integrals, mean slopes, sinuosity, 
%       Euclidean distance between the headwater and channel head, and 
%       along-basin distance from the headwater to the largest width.
%   tot_db_grids: Structure of grids for overall basin
%       lengths, heights, widths, orientations (phi), flow lengths, 
%       hypsometry integrals, mean slopes, sinuosities, and drainage area.
%   cross_db_stats: Array of drainage network cross-basin statistics. Rows
%       are basin ids, drainage sampling x-coordinates, drainage sampling 
%       y-coordinates, drainage sampling elevations, cross-basin widths 
%       (perpendicular to broad basin orientation) cross-basin reliefs 
%       (perpendicular to broad basin orientation).

%% Setup
% Get Drainage Basins
[DBg,x,y] = GRIDobj2mat(DB);
[X,Y] = meshgrid(x,y);
DBg = double(DBg);
DBg(isnan(Z)) = NaN;
dbi = unique(DBg(:));
dbi(isnan(dbi)) = [];

% Get flow directions
[Dg,~,~] = GRIDobj2mat(D);
Dg = double(Dg);
Dg(isnan(Z)) = NaN;

% Get slopes
[SlpG,~,~] = GRIDobj2mat(Slp);

% Get Drainage Areas
[Ag,~,~] = GRIDobj2mat(A);
Ag = Ag.*A.cellsize^2;

% Initialize values

tot_db_stats = zeros(length(dbi),12);
tot_db_grids = [];
cross_db_stats = [];

tmpWidths = DB;
tmpLengths = DB;
tmpHeights = DB;
tmpOrients = DB;
tmpHyps = DB;
tmpSlopes = DB;
tmpSinuosity = DB;
tmpFlowLengths = DB;
tmpDrainageArea = DB;
tmpEucDist = DB;

Xv = double(X(:));
Yv = double(Y(:));
Zv = double(Z(:));

%% Collect Stats on Each Basin
for i = 1:length(dbi)
    bL = DBg==dbi(i);
    
    % Get basin boundary points
    bb = bwboundaries(bL==1);
    boundXYZ = [];
    for j = 1:size(bb{1},1)
        boundXYZ = [boundXYZ;x(bb{1}(j,2)),y(bb{1}(j,1)),Z(bb{1}(j,1),bb{1}(j,2))];
    end
    boundXYZ = double(boundXYZ);
    
    % Get outlet and most-distant channel head.
    tmpD = Dg;
    tmpD(bL==0) = NaN;
    [maxDi,maxDj] = find(tmpD==max(tmpD(:)),1);
    
    % Extract channel from channel head to outlet.
    maxDLin = sub2ind(size(tmpD),maxDi,maxDj);
    [channelZi,~,channelX,channelY] = flowpathextract(FD,maxDLin);
    
    [inp,onp] = inpolygon(channelX,channelY,boundXYZ(:,1),boundXYZ(:,2));
    pp = inp+onp;
    channelX(pp==0) = [];
    channelY(pp==0) = [];
    channelZi(pp==0) = [];

    if isempty(channelX)
        channelX = 0;
        channelY = 0;
        channelZi = 1;
    end
    cHeadXY = [channelX(1),channelY(1)];
    cOutXY = [channelX(end),channelY(end)];
    channelZ = Z(channelZi);

    totEucDist = sqrt((cHeadXY(1)-cOutXY(1)).^2+(cHeadXY(2)-cOutXY(2)).^2);

    % Get Flow Length
    if length(channelX) > 1
        pp = [channelX,channelY];
        totFLength = diff(pp,1);
        totFLength = sum(sqrt(sum(totFLength.*totFLength,2)));
    else
        totFLength = NaN;
    end

    % Get basin hypsometry values
    tmpZ = Z;
    tmpZ(bL==0) = NaN;
    [bas_hyps_areas,bas_hyps_vals] = DrainageVolc_HypsometryValue(tmpZ,.01,1);
    totHyps = trapz(bas_hyps_areas,bas_hyps_vals);

    % Get local slope values
    tmpS = SlpG;
    tmpS(bL==0) = NaN;
    totSlopes = nanmean(tmpS(:));

    % Get upstream drainage area
    totDArea = max(Ag(bL==1));
    
    t1 = cHeadXY(1)==cOutXY(1);
    t2 = cHeadXY(2)==cOutXY(2);
    if t1*t2~=1
        %% Tot DB Stats & Grids
        % Get all basin points.
        dXs = X(bL);
        dYs = Y(bL);
        dZs = Z(bL);

        % Calcualte total basin height.
        totBHeight = range(dZs);

        % Calculate orientation by shifting outlet points relative to head
        % points and getting azimuth relative to the origin.
        dXs = cOutXY(1)-cHeadXY(1);
        dYs = cOutXY(2)-cHeadXY(2);
        sl = dYs/dXs;

        if dXs > 0
            if dYs > 0
                phi = 90-atand(abs(sl));
            elseif dYs < 0
                phi = 90+atand(abs(sl));
            else
                phi = 90;
            end
        elseif dXs < 0
            if dYs > 0
                phi = 270 + atand(abs(sl));
            elseif dYs < 0
                phi = 270 - atand(abs(sl));
            else
                phi = 270;
            end
        else
            if dYs > 0
                phi = 0;
            elseif dYs < 0
                phi = 180;
            else
                phi = NaN;
            end
        end

        totBOrient = phi;

        % Rotate Data
        mathPhi = 450-phi;
        if phi>=360
            phi = phi-360;
        end

        %% Get cross-basin stastics
        % Rotate all boundary points to channel head
        dXs = boundXYZ(:,1)-channelX(1);
        dYs = boundXYZ(:,2)-channelY(1);

        rotMat = [cosd((mathPhi)) -sind((mathPhi));sind((mathPhi)) cosd((mathPhi))];
        tmpYX = rotMat*[dYs,dXs]';
        nBoundXYs = double([tmpYX(2,:)',tmpYX(1,:)']);    

        % Rotate all channel points to channel head
        dXs = channelX-channelX(1);
        dYs = channelY-channelY(1);

        rotMat = [cosd((mathPhi)) -sind((mathPhi));sind((mathPhi)) cosd((mathPhi))];
        tmpYX = rotMat*[dYs,dXs]';
        nChannelXYs = [tmpYX(2,:)',tmpYX(1,:)'];
        
        % Rotate all DEM points to channel head
        gridYX = rotMat*[Yv-channelY(1),Xv-channelX(1)]';
        nGridXYs = [gridYX(2,:)',gridYX(1,:)'];
        gridInterp = scatteredInterpolant(nGridXYs(:,1),nGridXYs(:,2),Zv);

        % Run through channel points, get widths
        crossWidths = zeros(size(nChannelXYs(:,1)))*NaN;
        crossHeights = crossWidths;
        cross_midXY = zeros(size(nChannelXYs));
        if usePar
            parfor j = 1:length(crossWidths)
                xyDists = [[1:length(nBoundXYs(:,1))]',abs(nBoundXYs(:,1)-nChannelXYs(j,1)),nBoundXYs(:,2)-nChannelXYs(j,2)];
                xyDUpper = xyDists(xyDists(:,3)>=0,:);
                xyDLower = xyDists(xyDists(:,3)<0,:);
                
                if isempty(xyDUpper) || isempty(xyDLower)
                    continue;
                end
                
                upperI = find(abs(xyDUpper(:,2)) == min(abs(xyDUpper(:,2))),1);
                lowerI = find(abs(xyDLower(:,2)) == min(abs(xyDLower(:,2))),1);
                
                crossWidths(j) = abs(xyDUpper(upperI,3)) + abs(xyDLower(lowerI,3));
                cross_midXY(j,:) = [nChannelXYs(j,1),mean([nBoundXYs(xyDUpper(upperI,1),2);nBoundXYs(xyDLower(lowerI,1),2)],1)];
                
                rotBoundPoint1 = nBoundXYs(xyDUpper(upperI,1),:);
                rotBoundPoint2 = nBoundXYs(xyDLower(lowerI,1),:);
                
                ddx = rotBoundPoint1(1)-rotBoundPoint2(1);
                ddy = rotBoundPoint1(2)-rotBoundPoint2(2);
                slp = ddy/ddx;
                hyp = sqrt(ddx^2+ddy^2);
                newHyps = linspace(0,hyp,30);
                
                if ddx ~=0
                    interpX = rotBoundPoint2(1)+sign(ddx)*newHyps*cos(atan(abs(slp)));
                    interpY = rotBoundPoint2(2)+sign(ddy)*newHyps*sin(atan(abs(slp)));
                else
                    interpY = rotBoundPoint2(2)+newHyps;
                    interpX = ones(size(interpY))*rotBoundPoint2(1);
                end
                
                interpZ = gridInterp(interpX,interpY);
                crossHeights(j) = range(interpZ);
            end
        else
            for j = 1:length(crossWidths)
                xyDists = [[1:length(nBoundXYs(:,1))]',abs(nBoundXYs(:,1)-nChannelXYs(j,1)),nBoundXYs(:,2)-nChannelXYs(j,2)];
                xyDUpper = xyDists(xyDists(:,3)>=0,:);
                xyDLower = xyDists(xyDists(:,3)<0,:);
                
                if isempty(xyDUpper) || isempty(xyDLower)
                    continue;
                end
                
                upperI = find(abs(xyDUpper(:,2)) == min(abs(xyDUpper(:,2))),1);
                lowerI = find(abs(xyDLower(:,2)) == min(abs(xyDLower(:,2))),1);
                
                crossWidths(j) = abs(xyDUpper(upperI,3)) + abs(xyDLower(lowerI,3));
                cross_midXY(j,:) = [nChannelXYs(j,1),mean([nBoundXYs(xyDUpper(upperI,1),2);nBoundXYs(xyDLower(lowerI,1),2)],1)];
                
                rotBoundPoint1 = nBoundXYs(xyDUpper(upperI,1),:);
                rotBoundPoint2 = nBoundXYs(xyDLower(lowerI,1),:);
                
                ddx = rotBoundPoint1(1)-rotBoundPoint2(1);
                ddy = rotBoundPoint1(2)-rotBoundPoint2(2);
                slp = ddy/ddx;
                hyp = sqrt(ddx^2+ddy^2);
                newHyps = linspace(0,hyp,30);
                
                if ddx ~=0
                    interpX = rotBoundPoint2(1)+sign(ddx)*newHyps*cos(atan(abs(slp)));
                    interpY = rotBoundPoint2(2)+sign(ddy)*newHyps*sin(atan(abs(slp)));
                else
                    interpY = rotBoundPoint2(2)+newHyps;
                    interpX = ones(size(interpY))*rotBoundPoint2(1);
                end
                
                interpZ = gridInterp(interpX,interpY);
                crossHeights(j) = range(interpZ);
            end
        end

        % Smooth Mid-Basin Data
        cross_midXY = sortrows(cross_midXY,1);
        ddx = mean(diff(cross_midXY(:,1)));
        pixSize = ceil(smoothBasinPointWavelength/ddx);

        try
            cross_midXY(:,2) = smoothdata(cross_midXY(:,2),'movmean',pixSize);

            totBWidth = max(crossWidths);
            BWidthLoc = find(crossWidths==totBWidth,1);
            totBLength = diff(cross_midXY,1);
            BWidthDistance = sum(sqrt(sum(totBLength(1:BWidthLoc,:).^2,2)));
            totBLength = sum(sqrt(sum(totBLength.^2,2)));
            totSinuosity = totFLength./totBLength;
        catch er
            totBWidth = NaN;
            totBLength = NaN;
            totSinuosity = NaN;
            BWidthDistance = NaN;
        end
        
    else
        disp(sprintf('      Basin %d too small for analysis',dbi(i)))
        totBLength = NaN;
        totBWidth = NaN;
        totBOrient = NaN;
        totBHeight = NaN;
        totSinuosity = NaN;
        BWidthDistance = NaN;
        
        channelX = NaN;
        channelY = NaN;
        channelZ = NaN;
        crossWidths = NaN;
        crossHeights = NaN;
    end

    %% Fill stats matrix
    tot_db_stats(i,:) = [dbi(i),totDArea,totFLength,totBLength,totBWidth,totBHeight,totBOrient,totHyps,totSlopes,totSinuosity,totEucDist,BWidthDistance];

    % Fill DB grids
    tmpWidths.Z(DBg==dbi(i)) = totBWidth;
    tmpLengths.Z(DBg==dbi(i)) = totBLength;
    tmpHeights.Z(DBg==dbi(i)) = totBHeight;
    tmpOrients.Z(DBg==dbi(i)) = totBOrient;
    tmpFlowLengths.Z(DBg==dbi(i)) = totFLength;
    tmpHyps.Z(DBg==dbi(i)) = totHyps;
    tmpSlopes.Z(DBg==dbi(i)) = totSlopes;
    tmpSinuosity.Z(DBg==dbi(i)) = totSinuosity;
    tmpDrainageArea.Z(DBg==dbi(i)) = totDArea;
    tmpEucDist.Z(DBg==dbi(i)) = totEucDist;
    
    % Fill cross values
    cross_db_stats = [cross_db_stats;[ones(size(channelX))*i,...
        channelX,channelY,channelZ,crossWidths,crossHeights]];
end
% Put grids in structure
tot_db_grids.BasinWidths = tmpWidths;
tot_db_grids.BasinHeights = tmpHeights;
tot_db_grids.BasinLengths = tmpLengths;
tot_db_grids.BasinOrients = tmpOrients;
tot_db_grids.FlowLengths = tmpFlowLengths;
tot_db_grids.BasinHyps = tmpHyps;
tot_db_grids.BasinSlopes = tmpSlopes;
tot_db_grids.FlowSinuosity = tmpSinuosity;
tot_db_grids.DrainageArea = tmpDrainageArea;
tot_db_grids.BasinEuclideanLength = tmpEucDist;
end
