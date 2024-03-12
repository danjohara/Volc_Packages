function [distBins,BCounts_1,BCounts_2,Basin_Distance_Stats,nBCounts_1,nBCounts_2,BID_1,BID_2,NormR_1,NormR_2,nABCounts_1,nABCounts_2] = Calculate_BasinNum_ByDistance(mets,stepSize,volcName,makePlot,cutByThreshold)
    addpath('C:\Users\Daniel\Desktop\Research\Belgium\DEM_Analysis\Scripting_Packages\General_Scripts\')
    fs = 12;

    %% Increase boundary points
    bxy = mets.GeneralParams.boundaryXY;
    newBXY = [];
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

    %% Get radial distances
    cutZ = mets.GeneralParams.cutZ;
    cutX = mets.GeneralParams.cutX;
    cutY = mets.GeneralParams.cutY;

    [cutX,cutY] = meshgrid(cutX,cutY);

    [i,j] = find(cutZ==max(cutZ(:)),1);
    cutR = sqrt((cutX-cutX(i,j)).^2 + (cutY-cutY(i,j)).^2);

    cutR(isnan(cutZ)) = NaN;

    Br = sqrt((bxy(:,1)-cutX(i,j)).^2 + (bxy(:,2)-cutY(i,j)).^2);

    %% Calculate grid azimuths
    cutR_Norm1 = cutR./max(cutR(:));
    cutR_Norm2 = cutR_Norm1*0;

    dXs = cutX-cutX(i,j);
    dYs = cutY-cutY(i,j);
    sl = dYs./dXs;

    az = zeros(size(cutX));
    az((dXs>0).*(dYs>0) == 1) = 90-atand(abs(sl((dXs>0).*(dYs>0) == 1)));
    az((dXs>0).*(dYs<0) == 1) = 90+atand(abs(sl((dXs>0).*(dYs<0) == 1)));
    az((dXs>0).*(dYs==0) == 1) = 90;
    az((dXs<0).*(dYs>0) == 1) = 270 + atand(abs(sl((dXs<0).*(dYs>0) == 1)));
    az((dXs<0).*(dYs<0) == 1) = 270 - atand(abs(sl((dXs<0).*(dYs<0) == 1)));
    az((dXs<0).*(dYs==0) == 1) = 270;
    az((dXs==0).*(dYs>0) == 1) = 0;
    az((dXs==0).*(dYs<0) == 1) = 180;
    
    az(isnan(cutZ)) = NaN;
    az_r = [az(:),cutR(:)];
%     az_r(sum(isnan(az_r),2)>0,:) = [];

    %% Calculate boundary azimuths
    dXs = bxy(:,1)-cutX(i,j);
    dYs = bxy(:,2)-cutY(i,j);
    sl = dYs./dXs;

    Baz = zeros(size(bxy(:,1)));
    Baz((dXs>0).*(dYs>0) == 1) = 90-atand(abs(sl((dXs>0).*(dYs>0) == 1)));
    Baz((dXs>0).*(dYs<0) == 1) = 90+atand(abs(sl((dXs>0).*(dYs<0) == 1)));
    Baz((dXs>0).*(dYs==0) == 1) = 90;
    Baz((dXs<0).*(dYs>0) == 1) = 270 + atand(abs(sl((dXs<0).*(dYs>0) == 1)));
    Baz((dXs<0).*(dYs<0) == 1) = 270 - atand(abs(sl((dXs<0).*(dYs<0) == 1)));
    Baz((dXs<0).*(dYs==0) == 1) = 270;
    Baz((dXs==0).*(dYs>0) == 1) = 0;
    Baz((dXs==0).*(dYs<0) == 1) = 180;

    Baz_r = [Baz,Br];

    da = [0;diff(Baz_r(:,1))];
    ii = find(abs(da)>300,1);
    Baz_r = circshift(Baz_r,-ii+1,1);
    Baz_r = unique(Baz_r,'stable','rows');
    
    %% Interpolate boundary to points
    newBaz_r = [];
    azStep = 1;
    az = min(Baz_r(:,1)):azStep:max(Baz_r(:,1));

    for i = 1:length(az)
        t1 = Baz_r(:,1) >= az(i)-azStep/2;
        t2 = Baz_r(:,1) <= az(i)+azStep/2;
        t3 = t1.*t2;
        newBaz_r = [newBaz_r;az(i),max(Baz_r(t3==1,2))];
    end

%     Baz_r = unique(Baz_r,'rows');
%     maxRs = interp1(Baz_r(:,1),Baz_r(:,2),unique(az_r(:,1)),'nearest','extrap');
%     maxRs = reshape(maxRs,size(cutX));

    maxRs = interp1(newBaz_r(:,1),newBaz_r(:,2),az_r(:,1),'linear','extrap');
    cutR_Norm2 = cutR./reshape(maxRs,size(cutX));
    cutR_Norm2 = cutR_Norm2./max(cutR_Norm2(:));

    %% Perform counts
    distBins = .01:.01:1;
    spacingDistBins = .1:.1:.9;
    [DBg,~,~] = GRIDobj2mat(mets.DrainageParams.DB);
    [Ag,~,~] = GRIDobj2mat(mets.DrainageParams.Basin_Statistics_Grids.DrainageArea);
    if cutByThreshold == 1
        DBg(Ag<mets.ChannelParams.ChannelThreshold) = NaN;
    end
    BCounts_1 = zeros(size(distBins));
    BCounts_2 = BCounts_1;
    nBCounts_1 = BCounts_1;
    nBCounts_2 = BCounts_1;
    nABCounts_1 = BCounts_1;
    nABCounts_2 = BCounts_1;
    BID_1 = {};
    BID_2 = BID_1;

    for i = 1:length(distBins)
        tmp = DBg;
        tmp(cutR_Norm1>distBins(i)) = NaN;
        uniDB = unique(tmp(:));
        uniDB(isnan(uniDB)) = [];
        BCounts_1(i) = length(uniDB);
        BID_1{i} = uniDB;
        cc = contourc(cutX(1,:),cutY(:,1),cutR_Norm1,[1,1]*distBins(i));
        cc = Convert_Contours(cc,1);
        cumDist = 0;
        for j = 2:size(cc,1)
            cumDist = cumDist + sqrt((cc(j,1)-cc(j-1,1))^2 + (cc(j,2)-cc(j-1,2))^2);
        end
        cumDist = cumDist + sqrt((cc(1,1)-cc(end,1))^2 + (cc(1,2)-cc(end,2))^2);
        nBCounts_1(i) = length(uniDB)/cumDist;
        nABCounts_1(i) = length(uniDB)/polyarea(cc(:,1),cc(:,2));

%         if distBins(i) == .3
%             disp('here')
%         end
        tmp = DBg;
        tmp(cutR_Norm2>distBins(i)) = NaN;
        uniDB = unique(tmp(:));
        uniDB(isnan(uniDB)) = [];
        BCounts_2(i) = length(uniDB);
        BID_2{i} = uniDB;
        cc = contourc(cutX(1,:),cutY(:,1),cutR_Norm2,[1,1]*distBins(i));
        cc = Convert_Contours(cc,1);
        cumDist = 0;
        for j = 2:size(cc,1)
            cumDist = cumDist + sqrt((cc(j,1)-cc(j-1,1))^2 + (cc(j,2)-cc(j-1,2))^2);
        end
        cumDist = cumDist + sqrt((cc(1,1)-cc(end,1))^2 + (cc(1,2)-cc(end,2))^2);
        nBCounts_2(i) = length(uniDB)/cumDist;
        nABCounts_2(i) = length(uniDB)/polyarea(cc(:,1),cc(:,2));
    end

    boundCumDists = 0;
    for i = 2:size(bxy,1)
        tmpDist = sqrt((bxy(i,1)-bxy(i-1,1))^2 + (bxy(i,2)-bxy(i-1,2))^2);
        boundCumDists = [boundCumDists;boundCumDists(end)+tmpDist];
    end

    All_Basin_Dists1 = {};
    Average_Basin_Dists1 = [];
    STD_Basin_Dists1 = [];

    All_Basin_Dists2 = {};
    Average_Basin_Dists2 = [];
    STD_Basin_Dists2 = [];

    for i = 1:length(spacingDistBins)
        % Method 1
        tmpDB = DBg;
        tmpDB(cutR_Norm1>spacingDistBins(i)) = NaN;
        uniDB = unique(tmpDB);
        uniDB(isnan(uniDB)) = [];

        topBasinOutletPos = [];
        tmpDB = DBg;
        for j = 1:length(uniDB)
            tmpA = Ag;
            tmpA(tmpDB~=uniDB(j)) = NaN;
            [ii,jj] = find(tmpA==max(tmpA(:)),1);
            tmpDists = pdist2([cutX(ii,jj),cutY(ii,jj)],bxy);
            [~,bI] = find(tmpDists == min(tmpDists),1);
            topBasinOutletPos = [topBasinOutletPos;bI];
        end
        topBasinOutletPos = sortrows(topBasinOutletPos);
        basinBoundDists = [];
        try
        for j = 2:length(topBasinOutletPos)
            basinBoundDists = [basinBoundDists;
                boundCumDists(topBasinOutletPos(j))-boundCumDists(topBasinOutletPos(j-1))];
        end
        catch er
            disp('here')
        end
        basinBoundDists = [basinBoundDists;(boundCumDists(end)-boundCumDists(topBasinOutletPos(end)))+boundCumDists(topBasinOutletPos(1))];
    
        All_Basin_Dists1{i} = basinBoundDists;
        Average_Basin_Dists1(i) = mean(basinBoundDists);
        STD_Basin_Dists1(i) = std(basinBoundDists);

        % Method 2
        tmpDB = DBg;
        tmpDB(cutR_Norm2>spacingDistBins(i)) = NaN;
        uniDB = unique(tmpDB);
        uniDB(isnan(uniDB)) = [];

        topBasinOutletPos = [];
        tmpDB = DBg;
        for j = 1:length(uniDB)
            tmpA = Ag;
            tmpA(tmpDB~=uniDB(j)) = NaN;
            [ii,jj] = find(tmpA==max(tmpA(:)),1);
            tmpDists = pdist2([cutX(ii,jj),cutY(ii,jj)],bxy);
            [~,bI] = find(tmpDists == min(tmpDists),1);
            topBasinOutletPos = [topBasinOutletPos;bI];
        end
        topBasinOutletPos = sortrows(topBasinOutletPos);
        basinBoundDists = [];
        for j = 2:length(topBasinOutletPos)
            basinBoundDists = [basinBoundDists;
                boundCumDists(topBasinOutletPos(j))-boundCumDists(topBasinOutletPos(j-1))];
        end
        basinBoundDists = [basinBoundDists;(boundCumDists(end)-boundCumDists(topBasinOutletPos(end)))+boundCumDists(topBasinOutletPos(1))];
    
        All_Basin_Dists2{i} = basinBoundDists;
        Average_Basin_Dists2(i) = mean(basinBoundDists);
        STD_Basin_Dists2(i) = std(basinBoundDists);
    end

    Basin_Distance_Stats.Spacing_Distance_Bins = spacingDistBins;
    Basin_Distance_Stats.All_Basin_Distances_1 = All_Basin_Dists1;
    Basin_Distance_Stats.All_Basin_Distances_2 = All_Basin_Dists2;
    Basin_Distance_Stats.Average_Basin_Distances_1 = Average_Basin_Dists1;
    Basin_Distance_Stats.Average_Basin_Distances_2 = Average_Basin_Dists2;
    Basin_Distance_Stats.STD_Basin_Distances_1 = STD_Basin_Dists1;
    Basin_Distance_Stats.STD_Basin_Distances_2 = STD_Basin_Dists2;
    
    %% Make GRIDobj
    NormR_1 = GRIDobj(cutX,cutY,cutR_Norm1);
    NormR_2 = GRIDobj(cutX,cutY,cutR_Norm2);

    %% Make Plots
    if makePlot
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1)
        tmpDEM = GRIDobj(cutX,cutY,cutR_Norm1);
        topoDEM = GRIDobj(cutX,cutY,cutZ);
        imageschs(topoDEM,tmpDEM,'colormap',jet(255))
        hold on
        plot(bxy(:,1),bxy(:,2),'-k')
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        cb = colorbar;
        ylabel(cb,'Normalized Radial Distance')
        setAxes(cb,NaN)
        title({volcName,'(Method 1)'})
    
        subplot(1,2,2)
        tmpDEM = GRIDobj(cutX,cutY,cutR_Norm2);
        topoDEM = GRIDobj(cutX,cutY,cutZ);
        imageschs(topoDEM,tmpDEM,'colormap',jet(255))
        hold on
        plot(bxy(:,1),bxy(:,2),'-k')
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        cb = colorbar;
        ylabel(cb,'Normalized Radial Distance')
        setAxes(cb,NaN)
        title({volcName,'(Method 2)'})
    
    %     orient landscape
    %     set(gcf,'PaperType','tabloid')
    
        saveas(gcf,[volcName,'_Distances.fig'],'fig')
        saveas(gcf,[volcName,'_Distances.png'],'png')
    end
end

function setAxes(h,fs)
    set(h,'xcolor','k')
    set(h,'ycolor','k')
    if ~isnan(fs)
        set(h,'fontsize',fs)
    end
end