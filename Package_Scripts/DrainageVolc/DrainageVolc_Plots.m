function DrainageVolc_Plots(res)
%%
% Name: DrainageVolc_Plots
% Date: 06/30/2021 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to plot results of the DrainageVolc analysis.
%
% Input:
%   res: Result package from DrainageVolc.

%% Unpack Structure
GI = res.GeneralParams;
GP = res.GeographicParams;

DP = res.DrainageParams;
TP = res.TopoParams;
CP = res.ChannelParams;
DV = res.DivideParams;

useCM = viridis(255);

fs = 15;
lw = 3;
blw = 2;

visPlots = GI.inputs.visPlots;
tPre = GI.inputs.figTitlePrefix;
if ~isempty(tPre) && ~strcmp(tPre(end),' ')
    tPre = [tPre,' '];
end

[cutZ,cutx,cuty] = GRIDobj2mat(GP.DEM);
[cutX,cutY] = meshgrid(cutx,cuty);

kmS = CP.Channels.S;
kmS.x = kmS.x./1000;
kmS.y = kmS.y./1000;

kmTotChiS = CP.Chi.Total_Chi.Chi_S;
kmTotChiS.x = kmTotChiS.x./1000;
kmTotChiS.y = kmTotChiS.y./1000;

if isfield(DV.Divide_Topo,'DVD')
    kmDVD = DV.Divide_Topo.DVD;
    kmDVD.cellsize = kmDVD.cellsize/1000;
    kmDVD.distance = kmDVD.distance./1000;
    try
        kmDVD.refmat = kmDVD.refmat./1000;
    catch
        kmDVD.wf = kmDVD.wf./1000;
    end
else
    kmDVD = [];
end

%% Raw DEM
try
    figure('Name','Raw DEM','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)

    subplot(1,2,1)
    hold on
    imageschs(GP.DEM0,[],'colormap',demcmap(GP.DEM0.Z(:)),'tickstokm',true)
    allP = plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',blw);
    allPT = {'Edifice Boundary'};
    if ~isempty(GP.maskXY)
        for i = 1:length(GP.maskXY)
            p2 = fill(GP.maskXY{i}(:,1),GP.maskXY{i}(:,2),'w');
        end
        allP = [allP,p2];
        allPT = [allPT,{'Mask Region'}];
    end
    if ~isempty(GP.craterXY)
        for i = 1:length(GP.craterXY)
            p3 = fill(GP.craterXY{i}(:,1)./1000,GP.craterXY{i}(:,2)./1000,'r','facealpha',.5);
        end
        allP = [allP,p3];
        allPT = [allPT,{'Crater Region'}];
    end

    legend(allP,allPT);
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs);

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title([tPre,'Raw Elevations'])
    

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z(:)),'tickstokm',true)
    allP = plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',blw);
    allPT = {'Edifice Boundary'};
    if ~isempty(GP.maskXY)
        for i = 1:length(GP.maskXY)
            p2 = fill(GP.maskXY{i}(:,1),GP.maskXY{i}(:,2),'w');
        end
        allP = [allP,p2];
        allPT = [allPT,{'Mask Region'}];
    end
    if ~isempty(GP.craterXY)
        for i = 1:length(GP.craterXY)
            p3 = fill(GP.craterXY{i}(:,1)./1000,GP.craterXY{i}(:,2)./1000,'r','facealpha',.5);
        end
        allP = [allP,p3];
        allPT = [allPT,{'Crater Region'}];
    end

    legend(allP,allPT);
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs);

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title([tPre,'Clipped Elevations'])
    savePlot(GI,gcf,'Raw_Topography')
catch
    warning('Could not raw DEM')
end

%% Roughness
try
    if length(GI.inputs.roughnessWindows)>4
        warning('Number of roughness windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.roughnessWindows)]);

    figure('Name','Roughness Maps','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Roughness Maps
    for i = 1:maxSub
        gM = max(abs(TP.Roughness(i).Roughness_Grid.Z(:)));
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.Roughness(i).Roughness_Grid,'colormap',flipud(crameri('roma',255)),'caxis',[-1,1]*gM,'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Roughness (m)')
        setAxes(cb,0)
        title([tPre,sprintf('%.1f m Roughness',TP.Roughness(i).TrueWindowRes)])
    end

    savePlot(GI,gcf,'Roughness_Maps')
catch
    warning('Could not plot Roughness')
end

%% Windowed Slope Variance
try
    if length(GI.inputs.slopeVarianceWindows)>4
        warning('Number of slope variance windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.slopeVarianceWindows)]);

    figure('Name','Slope Variance Maps','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Slope Variance Maps
    for i = 1:maxSub
        gM = max(abs(TP.SlopeVariance_Windows(i).SlopeVariance_Grid.Z(:)));
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.SlopeVariance_Windows(i).SlopeVariance_Grid,'colormap',flipud(crameri('roma',255)),'caxis',[-1,1]*gM,'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Slope Variance')
        setAxes(cb,0)
        title({[tPre,sprintf('%.1f m',TP.SlopeVariance_Windows(i).TrueWindowRes)];'Slope Variance'})
    end

    savePlot(GI,gcf,'SlopeVariance_Windowed_Maps')
catch
    warning('Could not plot Windowed Slope Variance')
end

%% Drainage Area, Distance, and Basins
try
    figure('Name','Drainage Metrics','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Drainage Area Map
    subplot(2,2,1)
    hold on
    tmp = DP.Hydrology.A;
    tmp.Z = log10(tmp.Z*GI.inputs.dx^2);
    imageschs(GP.DEM,tmp,'colormap',flipud(crameri('tokyo',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,{'Log_{10} Cumulative';'Drainage Area (m^2)'})
    setAxes(cb,0)
    title([tPre,'Cumulative Drainage Area'])

    % Drainage Distance Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Hydrology.D,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs);

    cb = colorbar;
    ylabel(cb,'Flow Distance (m)')
    setAxes(cb,0)
    title([tPre,'Flow Distance'])

    % Drainage Basins Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.Drainage_Basins.DB,'colormap',crameri('hawaii',255),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Basin ID')
    setAxes(cb,0)
    title([tPre,'Drainage Basins'])

    % Drainage Basins Distribution
    subplot(2,2,4)
    hold on
    plot(DP.Hypsometry.DB_Hyps_numDB,DP.Hypsometry.DB_Hyps_Topo,'k','linewidth',2)
    plot(DP.Hypsometry.DB_Hyps_numDB_Channelized,DP.Hypsometry.DB_Hyps_Topo,'r','linewidth',2)
    box on
    set(gca,'xscale','log')
    setAxes(gca,fs)
    ylabel('Normalized Topography')
    xlabel('Number of Drainage Basins')
    axis tight
    axis square
    title([tPre,'Basin Topographic Distribution'])
    legend('All Basins','Channelized Basins','location','eastoutside')
    savePlot(GI,gcf,'Drainage_Metrics')
catch
    warning('Could not plot Drainage Values')
end

%% Hack's Law Plot - Basin Length
try
    tmpMap = DP.Hacks_Law.HackLawDeviation_BasinLength_Map;
    figure
    pcolor(tmpMap.Z);shading flat;
    bwrCM = colormap(bluewhitered(255));
    close

    HL = DP.Hacks_Law.HackLawDeviation_BasinLength;
    figure('Name','Hack''s Law Relationship (Basin Length)','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    colormap(bwrCM)
    subplot(1,2,1)
    hold on
    tmpRangeA = [min(HL(:,2)),max(HL(:,2))];
    tmp = scatter(HL(HL(:,5)==1,2),HL(HL(:,5)==1,3),ones(size(HL(HL(:,5)==1,2)))*200,HL(HL(:,5)==1,4),'o','filled','markeredgecolor','k');
    pp = [tmp];
    pL = {'Volcano Basins'};
    if sum(HL(:,5)==0) > 0
        tmp = plot(HL(HL(:,5)==0,2),HL(HL(:,5)==0,3),'.k');
        pp = [pp;tmp];
        pL = [pL;{'Excluded Basins'}];
        for i = 1:size(HL,1)
            if HL(i,5)==0
                tmpMap.Z(DP.Drainage_Basins.DB.Z==HL(i,1)) = NaN;
            end
        end
    end
    tmp = plot(tmpRangeA,DP.Hacks_Law.HackLawFit_BasinLength(1)*tmpRangeA.^DP.Hacks_Law.HackLawFit_BasinLength(2),'-k','linewidth',2);
    pp = [pp;tmp];
    pL = [pL;{sprintf('Best Fit (L = %.2f A^{%.2f})',DP.Hacks_Law.HackLawFit_BasinLength(1),DP.Hacks_Law.HackLawFit_BasinLength(2))}];
    xlim(tmpRangeA)

    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,fs)

    xlabel('Basin Drainage Area (m^2)')
    ylabel('Basin Length (m)')
    axis tight
    axis square
    box on
    if GI.inputs.limitHacksLaw
        yy = ylim;
        tmp = plot([1,1]*CP.Channels.ChannelThreshold,yy,'--r','linewidth',2);
        pp = [pp;tmp];
        pL = [pL;{'Channelization Threshold'}];
    end
    legend(pp,pL,'location','northwest')
    title([tPre,'Hack''s Law Relationship (Basin Length)'])

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,tmpMap,'colormap',bwrCM,'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
    plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title([tPre,sprintf('Basin Length Deviation from\nHack''s Law (Basin Length)')])
    savePlot(GI,gcf,'Hacks_Law_Relationship_BasinLength')
catch
    warning('Could not plot Basin Length Hack''s Law')
end

%% Hack's Law Plot - Basin Flow Length
try
    figure
    pcolor(DP.Hacks_Law.HackLawDeviation_FlowLength_Map.Z);shading flat;
    bwrCM = colormap(bluewhitered(255));
    close

    HL = DP.Hacks_Law.HackLawDeviation_FlowLength;
    figure('Name','Hack''s Law Relationship (Basin Flow Length)','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    colormap(bwrCM)
    subplot(1,2,1)
    hold on
    tmpRangeA = [min(HL(:,2)),max(HL(:,2))];
    tmp = scatter(HL(HL(:,5)==1,2),HL(HL(:,5)==1,3),ones(size(HL(HL(:,5)==1,2)))*200,HL(HL(:,5)==1,4),'o','filled','markeredgecolor','k');
    pp = [tmp];
    pL = {'Volcano Basins'};
    if sum(HL(:,5)==0) > 0
        tmp = plot(HL(HL(:,5)==0,2),HL(HL(:,5)==0,3),'.k');
        pp = [pp;tmp];
        pL = [pL;{'Excluded Basins'}];
    end
    tmp = plot(tmpRangeA,DP.Hacks_Law.HackLawFit_FlowLength(1)*tmpRangeA.^DP.Hacks_Law.HackLawFit_FlowLength(2),'-k','linewidth',2);
    pp = [pp;tmp];
    pL = [pL;{sprintf('Best Fit (L = %.2f A^{%.2f})',DP.Hacks_Law.HackLawFit_FlowLength(1),DP.Hacks_Law.HackLawFit_FlowLength(2))}];
    xlim(tmpRangeA)

    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,fs)

    xlabel('Basin Drainage Area (m^2)')
    ylabel('Max Basin Flow Length (m)')
    axis tight
    axis square
    box on
    if GI.inputs.limitHacksLaw
        yy = ylim;
        tmp = plot([1,1]*CP.Channels.ChannelThreshold,yy,'--r','linewidth',2);
        pp = [pp;tmp];
        pL = [pL;{'Channelization Threshold'}];
    end
    legend(pp,pL,'location','northwest')
    title([tPre,'Hack''s Law Relationship (Basin Flow Length)'])

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,DP.Hacks_Law.HackLawDeviation_FlowLength_Map,'colormap',bwrCM,'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
    plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title([tPre,sprintf('Basin Length Deviation from\nHack''s Law (Flow Length)')])
    savePlot(GI,gcf,'Hacks_Law_Relationship_FlowLength')
catch
    warning('Could not plot Flow Length Hack''s Law')
end

%% Drainage Area - Slope Plots
try
    figure('Name','Slope-Area Plots','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    subplot(1,2,1)
    plot(DP.Top_Drainage_Basins.TopN_All_Area_Slope(:,1),tand(DP.Top_Drainage_Basins.TopN_All_Area_Slope(:,2)),'.k')

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    setAxes(gca,15)
    yy = ylim();
    xx = xlim();
    xlabel('Drainage Area (m^2)')
    ylabel('Slope')
    axis square
    box on
    if GI.inputs.basinTopN >= 1
        title([tPre,sprintf('All Slope-Area of\nLargest %d Basins',length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    elseif GI.inputs.basinTopN > 0
        title([tPre,sprintf('All Slope-Area of\nLargest %d%% Basins\n(%d Basins)',GI.inputs.basinTopN*100,length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    elseif GI.inputs.basinTopN < 0 && GI.inputs.basinTopN > -1
        title([tPre,sprintf('All Slope-Area of Basins\nin Upper %d%% of Topography\n(%d Basins)',abs(GI.inputs.basinTopN)*100,length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    end

    subplot(1,2,2)
    hold on
    p1 = plot(DP.Top_Drainage_Basins.TopN_All_Area_Slope(:,1),tand(DP.Top_Drainage_Basins.TopN_All_Area_Slope(:,2)),'.','markeredgecolor',[.8 .8 .8]);
    p2 = plot(DP.Top_Drainage_Basins.TopN_Flow_Area_Slope(:,1),tand(DP.Top_Drainage_Basins.TopN_Flow_Area_Slope(:,2)),'.k');

    if ~isempty(DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs)
        as = DP.Top_Drainage_Basins.TopN_Flow_Area_Slope;
        as(as(:,2)==0,:) = [];
        as(as(:,1)<DP.Top_Drainage_Basins.TopN_TransitionThreshold_A,:) = [];
        AS1 = as(as(:,1)<=DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1),:);
        AS2 = as(as(:,1)>DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1),:);

        mdl1 = fitlm(log10(AS1(:,1)),log10(tand(AS1(:,2))));
        mdl2 = fitlm(log10(AS2(:,1)),log10(tand(AS2(:,2))));

        aa1 = [DP.Top_Drainage_Basins.TopN_TransitionThreshold_A,DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1)];
        aa2 = [DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1),max(as(:,1))];
        p3 = plot(aa1,10^mdl1.Coefficients.Estimate(1)*aa1.^mdl1.Coefficients.Estimate(2),'-r','linewidth',3);
        p4 = plot(aa2,10^mdl2.Coefficients.Estimate(1)*aa2.^mdl2.Coefficients.Estimate(2),'-b','linewidth',3);

        ss1 = tand(as(:,2));

        p5 = plot([DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1),DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1)],[min(ss1),max(ss1)],'--k','linewidth',3);
        p6 = plot([1,1]*DP.Top_Drainage_Basins.TopN_TransitionThreshold_A,[min(ss1),max(ss1)],'--r','linewidth',3);

        ps = [p1,p2,p6,p3,p4,p5];
        pt = {'All Slope-Area','Flowpath Slope-Area','Transition Zone Start',...
        sprintf('Regression 1 (r^2 = %.2f)',DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(2)),...
        sprintf('Regression 2 (r^2 = %.2f)(M/N = %.2f)',DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(3),DP.Top_Drainage_Basins.TopN_MN),...
        sprintf('Area Threhold = %.2f km^2 (r'' = %.2f)',DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(1)./1e6,DP.Top_Drainage_Basins.TopN_AreaThreshold_Rs(4))};

    else
        ps = [p1,p2];
        pt = {'All A-S','Flowpath A-S'};
    end

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    setAxes(gca,15)
    xlim(xx)
    ylim(yy)
    xlabel('Drainage Area (m^2)')
    ylabel('Slope')
    axis square
    box on
    legend(ps,pt,'location','southwest')
    if GI.inputs.basinTopN >= 1
        title([tPre,sprintf('Flowpath Slope-Area of\nLargest %d Basins',length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    elseif GI.inputs.basinTopN > 0
        title([tPre,sprintf('Flowpath Slope-Area of\nLargest %d%% Basins\n(%d Basins)',GI.inputs.basinTopN*100,length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    elseif GI.inputs.basinTopN < 0 && GI.inputs.basinTopN > -1
        title([tPre,sprintf('Flowpath Slope-Area of Basins\nin Upper %d%% of Topography\n(%d Basins)',abs(GI.inputs.basinTopN)*100,length(CP.Channel_Concavity.Concavity_BasinIDs))]);
    end
    savePlot(GI,gcf,'Slope-Area')
catch
    warning('Could not plot Slope-Area')
end

%% Basin Statistics 1
try
    figure('Name','Basin Total Statistics 1','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Basin Height Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinHeights,'colormap',flipud(crameri('acton',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title([tPre,'Basin Relief'])

    % Basin Length Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinLengths,'colormap',flipud(crameri('turku',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Length (m)')
    setAxes(cb,0)
    title([tPre,'Basin Length'])

    % Basin Width Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinWidths,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title([tPre,'Maximum Basin Width'])

    % Basin Ellipticity Map
    tmp = DP.Basin_Statistics.Basin_Statistics_Grids.BasinWidths;
    tmp.Z = tmp.Z./DP.Basin_Statistics.Basin_Statistics_Grids.BasinLengths.Z;
    tmp.Z(tmp.Z>1) = NaN;
    subplot(2,2,4)
    hold on
    imageschs(GP.DEM,tmp,'colormap',flipud(crameri('lajolla',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Ellipticity')
    setAxes(cb,0)
    title([tPre,'Basin Ellipticity'])

    savePlot(GI,gcf,'Basin_Total_Statistics_1')
catch
    warning('Could not plot Basin Statistics 1')
end

%% Basin Statistics 2
try
    figure('Name','Basin Total Statistics 2','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Basin Hypsometry
    subplot(2,2,1)
    hold on
    aa = DP.Basin_Statistics.Basin_Statistics_Grids.BasinHyps;
    aa.Z = abs(aa.Z);
    imageschs(GP.DEM,aa,'colormap',crameri('buda',255),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Hypsometry Integral')
    setAxes(cb,0)
    title([tPre,'Basin Hypsometry'])

    % Basin Mean Slope
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinSlopes,'colormap',flipud(crameri('imola',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    title([tPre,'Basin Mean Slopes'])

    % Drainage Area Map
    subplot(2,2,3)
    hold on
    tmp = DP.Basin_Statistics.Basin_Statistics_Grids.DrainageArea;
    tmp.Z = log10(tmp.Z);
    imageschs(GP.DEM,tmp,'colormap',flipud(crameri('tokyo',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Log_{10} Drainage Area (m^2)')
    setAxes(cb,0)
    title([tPre,'Basin Drainage Area'])

    % Basin Orientation Map
    subplot(2,2,4)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinOrients,'colormap',crameri('brocO',255),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Orientation (Azimuth)')
    setAxes(cb,0)
    title([tPre,'Basin Orientation'])

    savePlot(GI,gcf,'Basin_Total_Statistics_2')
catch
    warning('Could not plot Basin Statistics 2')
end
%% Basin Statistics 3
try
    figure('Name','Basin Total Statistics 3','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Basin Flow Length Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.FlowLengths,'colormap',flipud(crameri('lapaz',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Flow Length (m)')
    setAxes(cb,0)
    title([tPre,'Basin Flow Length'])

    % Basin Euclidean Length Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.BasinEuclideanLength,'colormap',flipud(crameri('lapaz',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Euclidean Length (m)')
    setAxes(cb,0)
    title([tPre,'Basin Euclidean Length'])
    cax = caxis;
    caxis([1,cax(2)])

    % Basin Sinuosity Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics.Basin_Statistics_Grids.FlowSinuosity,'colormap',crameri('bilbao',255),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Sinuosity')
    setAxes(cb,0)
    title([tPre,'Basin Flow Sinuosity'])
    cax = caxis;
    caxis([1,cax(2)])

    % Basin Order
    if ~isempty(CP.Channels.S)
        mo = nanmax(CP.Channels.Max_Stream_Order.Z(:));
        mSpacing = (mo-1)/mo;
        mStart = 1+mSpacing/2;
        subplot(2,2,4)
        hold on
        imageschs(GP.DEM,CP.Channels.Max_Stream_Order,'colormap',crameri('hawaii',mo),'tickstokm',true);
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Maximum Stream Order')
        caxis([1,mo])
        set(cb,'ytick',mStart:mSpacing:mo)
        set(cb,'yticklabel',1:mo)
    
        setAxes(cb,0)
        title([tPre,'Basin Strahler Stream Order'])
    end

    savePlot(GI,gcf,'Basin_Total_Statistics_3')
catch
    warning('Could not plot Basin Statistics 3')
end

%% Basin Sinuosity
try
    figure('Name','Basin Contour Sinuosity','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    tmpMeans = GP.DEM;
    tmpMeans.Z(:) = NaN;
    tmpMedians = tmpMeans;
    tmpMaxes = tmpMeans;
    tmpMins = tmpMeans;
    for i = 1:length(DP.Basin_Contour_Sinuosity_Stats.BasinIDs)
        tmpMeans.Z(DP.Drainage_Basins.DB.Z==DP.Basin_Contour_Sinuosity_Stats.BasinIDs(i)) = DP.Basin_Contour_Sinuosity_Stats.Means(i);
        tmpMedians.Z(DP.Drainage_Basins.DB.Z==DP.Basin_Contour_Sinuosity_Stats.BasinIDs(i)) = DP.Basin_Contour_Sinuosity_Stats.Medians(i);
        tmpMaxes.Z(DP.Drainage_Basins.DB.Z==DP.Basin_Contour_Sinuosity_Stats.BasinIDs(i)) = DP.Basin_Contour_Sinuosity_Stats.Maxes(i);
        tmpMins.Z(DP.Drainage_Basins.DB.Z==DP.Basin_Contour_Sinuosity_Stats.BasinIDs(i)) = DP.Basin_Contour_Sinuosity_Stats.Mins(i);
    end

    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,tmpMins,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    contour(cutX./1000,cutY./1000,cutZ,DP.Basin_Contour_Sinuosity_Stats.Contours,'-r','linewidth',1)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Contour Sinuosity')
    setAxes(cb,0)
    title([tPre,'Minimum Contour Sinuosity'])

    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,tmpMaxes,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    contour(cutX./1000,cutY./1000,cutZ,DP.Basin_Contour_Sinuosity_Stats.Contours,'-r','linewidth',1)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Contour Sinuosity')
    setAxes(cb,0)
    title([tPre,'Maximum Contour Sinuosity'])

    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,tmpMeans,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    contour(cutX./1000,cutY./1000,cutZ,DP.Basin_Contour_Sinuosity_Stats.Contours,'-r','linewidth',1)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Contour Sinuosity')
    setAxes(cb,0)
    title([tPre,'Mean Contour Sinuosity'])

    subplot(2,2,4)
    hold on
    imageschs(GP.DEM,tmpMedians,'colormap',flipud(crameri('davos',255)),'tickstokm',true)
    contour(cutX./1000,cutY./1000,cutZ,DP.Basin_Contour_Sinuosity_Stats.Contours,'-r','linewidth',1)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Contour Sinuosity')
    setAxes(cb,0)
    title([tPre,'Median Contour Sinuosity'])

    savePlot(GI,gcf,'Basin_Contour_Sinuosity')
catch
    warning('Could not plot Basin Contour Sinuosity')
end

%% Basin Roughness
try
    figure('Name','Mean Basin Roughness','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)

    if length(GI.inputs.roughnessWindows)>4
        warning('Number of roughness windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.roughnessWindows)]);

    for i = 1:maxSub
        tmpMeans = GP.DEM;
        tmpMeans.Z(:) = NaN;
        for j = 1:length(DP.Basin_Roughness_Stats.BasinIDs)
            tmpMeans.Z(DP.Drainage_Basins.DB.Z==DP.Basin_Roughness_Stats.BasinIDs(j)) = DP.Basin_Roughness_Stats.Means(j,i);
        end

        rMax = abs(tmpMeans.Z(:));
        rMax = nanmax(rMax);
        cmRange = [-1,1]*rMax;


        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,tmpMeans,'colormap',flipud(crameri('roma',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
        caxis(cmRange)

        cb = colorbar;
        ylabel(cb,'Mean Roughness')
        setAxes(cb,0)
        title([tPre,sprintf('Mean Basin Roughness\n(Window Size = %.2f km)',DP.Basin_Roughness_Stats.Windows(i)/1000)]);
    end

    savePlot(GI,gcf,'Basin_Mean_Roughness')
catch
    warning('Could not plot Basin Mean Roughness')
end

%% Basin Slope Variance
try
    figure('Name','Mean Basin Slope Variance','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)

    if length(GI.inputs.slopeVarianceWindows)>4
        warning('Number of slope variance windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.slopeVarianceWindows)]);

    for i = 1:maxSub
        tmpMeans = GP.DEM;
        tmpMeans.Z(:) = NaN;
        for j = 1:length(DP.Basin_SlopeVariance_Stats.BasinIDs)
            tmpMeans.Z(DP.Drainage_Basins.DB.Z==DP.Basin_SlopeVariance_Stats.BasinIDs(j)) = DP.Basin_SlopeVariance_Stats.Means(j,i);
        end

        rMax = abs(tmpMeans.Z(:));
        rMax = nanmax(rMax);
        cmRange = [-1,1]*rMax;


        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,tmpMeans,'colormap',flipud(crameri('roma',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
        caxis(cmRange)

        cb = colorbar;
        ylabel(cb,'Mean Slope Variance')
        setAxes(cb,0)
        title([tPre,sprintf('Mean Basin Slope Variance\n(Window Size = %.2f km)',DP.Basin_Roughness_Stats.Windows(i)/1000)]);
    end

    savePlot(GI,gcf,'Basin_Mean_SlopeVariance')
catch
    warning('Could not plot Basin Mean Slope Variance')
end

%% Cross-Basin Statistics
try
    figure('Name','Basin Cross Statistics','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    % Basin Widths
    subplot(2,2,1)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
    hold on
    plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5)
    scatter(DP.Basin_Statistics.Basin_Cross_Statistics(:,2)./1000,DP.Basin_Statistics.Basin_Cross_Statistics(:,3)./1000,[],DP.Basin_Statistics.Basin_Cross_Statistics(:,5),'o','filled')
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title([tPre,'Cross-Basin Widths'])
    colormap(gca,flipud(crameri('davos')))

    % Basin Relief
    subplot(2,2,2)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
    hold on
    plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5)
    scatter(DP.Basin_Statistics.Basin_Cross_Statistics(:,2)./1000,DP.Basin_Statistics.Basin_Cross_Statistics(:,3)./1000,[],DP.Basin_Statistics.Basin_Cross_Statistics(:,6),'o','filled')
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title([tPre,'Cross-Basin Relief'])
    colormap(gca,flipud(crameri('acton')))

    % Basin Incision
    subplot(2,2,[3,4])
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
    hold on
    plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5)
    scatter(DP.Basin_Statistics.Basin_Cross_Statistics(:,2)./1000,DP.Basin_Statistics.Basin_Cross_Statistics(:,3)./1000,[],log10(DP.Basin_Statistics.Basin_Cross_Statistics(:,7)),'o','filled')
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    caxis([-3,0])
    cb = colorbar;
    ylabel(cb,'Log_{10} Incision Index')
    setAxes(cb,0)
    title({[tPre,'Cross-Basin Incision Index']; '(Relief / Width)'})
    colormap(gca,flipud(crameri('lajolla')))

    savePlot(GI,gcf,'Cross_Basin_Statistics')
catch
    warning('Could not plot Cross-Basin Statistics')
end

%% Drainage Basins Per Contour
try
    figure('Name','Basins Per Contour','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    subplot(1,2,1)
    hold on
    bCCpCLA = DP.Drainage_Basins.Basin_Contour_ContourP_Count_Length_Area;
    bCCpCLA_C = DP.Drainage_Basins.Basin_Contour_ContourP_Count_Length_Area_Channelized;
    p1 = plot(bCCpCLA(:,3),bCCpCLA(:,1),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw);
    p2 = plot(bCCpCLA_C(:,3),bCCpCLA_C(:,1),'-or','markerfacecolor','r','markersize',10,'linewidth',blw);
    xlabel('Number of Basins')
    ylabel('Contour (m)')
    setAxes(gca,fs)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])
    yyaxis right
    plot(bCCpCLA(:,3),bCCpCLA(:,2),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw)
    setAxes(gca,fs)
    if DP.Top_Drainage_Basins.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Top_Drainage_Basins.Basin_TopN+1),'--r','linewidth',lw)
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    legend([p1,p2],{'All Basins','Channelized Basins'})
    title([tPre,'Number of Basins per Contour'])
    box on
    axis square

    subplot(1,2,2)
    hold on
    bCCpCLA = DP.Drainage_Basins.Basin_Contour_ContourP_Count_Length_Area;
    p1 = plot(bCCpCLA(:,3)./bCCpCLA(:,4)*1000,bCCpCLA(:,1),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw);
    p2 = plot(bCCpCLA_C(:,3)./bCCpCLA_C(:,4)*1000,bCCpCLA_C(:,1),'-or','markerfacecolor','r','markersize',10,'linewidth',blw);
    xlabel('Basins per Contour Length (km^{-1})')
    setAxes(gca,fs)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])

    yyaxis right
    plot(bCCpCLA(:,3)./bCCpCLA(:,4)*1000,bCCpCLA(:,2),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw)

    setAxes(gca,fs)
    title([tPre,'Number of Basins per Contour Length'])
    if DP.Top_Drainage_Basins.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Top_Drainage_Basins.Basin_TopN+1),'--r','linewidth',lw)
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    ylabel('Percent Relief')
    legend([p1,p2],{'All Basins','Channelized Basins'})
    box on
    axis square
    savePlot(GI,gcf,'Basin_per_Contour')
catch
    warning('Could not plot Basin-Contour Relationships')
end

%% Radial Basin Count
try
    [ii,jj] = find(cutZ == max(cutZ(:)),1);
    figure('Name','Radial Basin Count','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,DP.Radial_Analysis.Normalized_Radial_Distances,'colormap',(crameri('buda',255)),'tickstokm',true)
    plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
    plot(cutx(jj)/1000,cuty(ii)/1000,'ok','markerfacecolor','r','markersize',12)

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Normalized Distance')
    setAxes(cb,0)
    title([tPre,'Normalized Radial Distance'])

    subplot(1,2,2)
    hold on
    plot(DP.Radial_Analysis.Basin_Count(:,2),DP.Radial_Analysis.Basin_Count(:,3),'-k','linewidth',lw)
    plot(DP.Radial_Analysis.Basin_Count(:,2),DP.Radial_Analysis.Basin_Count(:,4),'-r','linewidth',lw)

    xlabel('Normalized Radial Distance')
    ylabel('Basin Count')
    setAxes(gca,fs)
    legend('All Basins','Channelized Basins','location','northwest')
    title([tPre,'Number of Radial Basins'])
    box on
    axis square
    savePlot(GI,gcf,'Basin_per_RadialDistance')
catch
    warning('Could not plot Radial basin analysis')
end

%% Channels
try
    if isempty(CP.Channels.S)
        warning('Cannot plot Channels; no channels exist')
    else
        figure('Name','Channels','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plot(kmS,'-k','linewidth',2)
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Elevation (m)')
        setAxes(cb,0)
        title([tPre,sprintf('Drainage Channels (Threshold = %.2f km^2)',CP.Channels.ChannelThreshold./1e6)])
    
        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false,'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plot(kmS,'-k','linewidth',2)
        scatter(CP.Channels.Knickpoint_ID_BID_XY_StreamDist_DZ(:,3)./1000,...
            CP.Channels.Knickpoint_ID_BID_XY_StreamDist_DZ(:,4)./1000,...
            ones(size(CP.Channels.Knickpoint_ID_BID_XY_StreamDist_DZ(:,3)))*60,...
            CP.Channels.Knickpoint_ID_BID_XY_StreamDist_DZ(:,6),'o','filled',...
            'markeredgecolor','k')
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Knickpoint Magnitude (m)')
        setAxes(cb,0)
        title([tPre,'Knickpoint Locations'])
        colormap(gca,flipud(crameri('nuuk')))
    
        savePlot(GI,gcf,'Channel')
    end
catch
    warning('Could not plot Channels')
end

%% Drainage Density Along Channels
try
    if isempty(CP.Channels.S)
        warning('Cannot plot Drainage Density Along Channels; no channels exist')
    else
        f = figure;
        cm = colormap(crameri('bilbao',255));
        cm = cm(1:end,:);
        close(f);
        
        figure('Name','Channel Drainage Density','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        % Channel Drainage Density Map
        hold on
        imageschs(GP.DEM,[],'colormap',[1 1 1]*.75,'colorbar',false,'tickstokm',true)
        plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plotc(kmS,CP.Drainage_Density.DD*1000,'linewidth',5); colormap(cm)
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Drainage Density (km^{-1})')
        setAxes(cb,0)
        box on
        title([tPre,'Drainage Density Along Channels'])
    
        savePlot(GI,gcf,'Drainage_Density_Channels')
    end
catch
    warning('Could not plot Drainage Density Along Channels')
end

%% Basin Drainage Density
try
    if isempty(CP.Channels.S)
        warning('Cannot plot Basin Drainage Density; no channels exist')
    else
        figure('Name','Basin Drainage Density','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        % Basin Drainage Density Map
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,CP.Drainage_Density.BasinDD*1000,'colormap',flipud(crameri('acton',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plot(kmS,'-k','linewidth',blw)
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Total Basin Drainage Density (km^{-1})')
        setAxes(cb,0)
        title([tPre,sprintf('Total Basin\nDrainage Density')])
    
        % Max Basin Drainage Density Map
        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,CP.Drainage_Density.MaxDD*1000,'colormap',flipud(crameri('oslo',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plot(kmS,'-k','linewidth',blw)
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Max Basin Drainage Density (km^{-1})')
        setAxes(cb,0)
        title([tPre,sprintf('Maximum Basin\nDrainage Density')])
        savePlot(GI,gcf,'Drainage_Density_Basins')
    end
catch
    warning('Could not plot Basin Drainage Density')
end

%% Channel Conformities
try
    if isempty(CP.Channels.S)
        warning('Cannot plot Channel Conformity; no channels exist')
    else
        figure('Name','Channel Conformity','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        subplot(2,2,1)
        hold on
        imageschs(GP.DEM,CP.Conformity.Filtered_Topography,'colormap',demcmap(GP.DEM.Z),'tickstokm',true)
        plot(kmS,'-k','linewidth',blw)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Elevation (m)')
        setAxes(cb,0)
        title([tPre,'Filtered Topography'])
    
        subplot(2,2,3)
        hold on
        imageschs(GP.DEM,CP.Conformity.Mean_Basin_Conformity_Map,'colormap',(crameri('acton',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        plot(kmS,'-k','linewidth',blw)
    
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
    
        cb = colorbar;
        ylabel(cb,'Conformity')
        setAxes(cb,0)
        title([tPre,sprintf('Mean Basin Conformity\n(Total Conformity = %.2f)',CP.Conformity.Mean_Total_Conformity)])
    
        subplot(2,2,[2,4])
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
    
        for i = 1:length(CP.Conformity.Segment_Information)
            xx = [CP.Conformity.Segment_Information(i).X,CP.Conformity.Segment_Information(i).X];
            yy = [CP.Conformity.Segment_Information(i).Y,CP.Conformity.Segment_Information(i).Y];
            zz = zeros(size(xx));
            cc = ones(size(xx))*CP.Conformity.Segment_Information(i).Conformity;
    
            surf(xx./1000,yy./1000,zz,cc,'edgecolor','interp','facecolor','none','linewidth',3)
            quiver(CP.Conformity.Segment_Information(i).LongWavelength_EndPoints(end,1)/1000,...
                CP.Conformity.Segment_Information(i).LongWavelength_EndPoints(end,2)/1000,...
                -diff(CP.Conformity.Segment_Information(i).LongWavelength_EndPoints(:,1))./1000,...
                -diff(CP.Conformity.Segment_Information(i).LongWavelength_EndPoints(:,2))./1000,'-k','linewidth',1)
        end
        colormap(gca,(crameri('acton')))
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
        cb = colorbar;
        ylabel(cb,'Conformity')
        setAxes(cb,0)
        title([tPre,'Segment Conformity'])
        box on
    
        savePlot(GI,gcf,'Basin_Conformity')
    end
catch
    warning('Could not plot Channel Conformity')
end
    
%% Channel Concavity
try
    if isempty(CP.Channels.S)
        warning('Cannot plot Channel Concavity; no channels exist')
    else
        figure('Name','Concavity','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
    
        if size(CP.Channel_Concavity.Concavity_Streams,1) == 0
            useBasinIs = [];
        else
            nonNaNIs = 1:size(CP.Channel_Concavity.Concavity_Stats,1);
            asLengths = zeros(size(nonNaNIs));
            t = ones(1,size(CP.Channel_Concavity.Concavity_Stats,1));
            for i = 1:size(CP.Channel_Concavity.Concavity_Stats,1)
                t(i) = isnan(CP.Channel_Concavity.Concavity_Stats(i).ks);
                asLengths(i) = length(CP.Channel_Concavity.Concavity_Stats(i).a);
            end
            nonNaNIs(t==1) = [];
            asLengths(t==1) = [];
            IL = sortrows([nonNaNIs;asLengths]',2,'descend');
    
            if length(nonNaNIs) <= 4
                useBasinIs = length(nonNaNIs);
            else
                % useBasinIs = randsample(nonNaNIs,4);
                useBasinIs = IL(1:4,1);
            end
        end
        useBasinIDs = [];
        % Concavity Plot
        for i = 1:length(useBasinIs)
            if i > 2
                useI = i+2;
            else
                useI = i;
            end
            subplot(2,4,useI)
            hold on
            cc = [CP.Channel_Concavity.Concavity_Stats(useBasinIs(i)).a,CP.Channel_Concavity.Concavity_Stats(useBasinIs(i)).g];
            cc(cc(:,2)<1e-3,:) = [];
            plot(cc(:,1),cc(:,2),'sk','markerfacecolor','r','markersize',10)
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            xx = xlim();
            yy = ylim();
            plot(sortrows(cc(:,1)),CP.Channel_Concavity.Concavity_Stats(i).ks(1)*sortrows(cc(:,1)).^CP.Channel_Concavity.Concavity_Stats(i).theta,'--k','linewidth',lw);
            xlim(xx);
            ylim(yy);
            setAxes(gca,fs)
            title([sprintf('Basin %d\n\\theta = %.2f',CP.Channel_Concavity.Concavity_BasinIDs(useBasinIs(i)),CP.Channel_Concavity.Concavity_Stats(useBasinIs(i)).theta)])
            box on
            if i > 2
                xlabel('Area (m^2)')
            end
            if i == 1 || i == 5
                ylabel('Slope')
            end
            axis square
            useBasinIDs = [useBasinIDs;CP.Channel_Concavity.Concavity_BasinIDs(useBasinIs(i))];
        end
    
        % Concavity Basin Map
        subplot(2,4,[3,4,7,8])
        hold on
        try
            imageschs(GP.DEM,CP.Channel_Concavity.Concavity_DEM,'colormap',flipud(crameri('devon',255)),'caxis',[0,2],'tickstokm',true)
        catch
            warning('No Concavity values available - decrease the drainage area threshold.');
        end
    
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',.5)
        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)
        title([tPre,'Drainage Basins'])
        [DBg,~,~] = GRIDobj2mat(DP.Drainage_Basins.DB);
        for i = 1:size(CP.Channel_Concavity.Concavity_Streams,1)
            t1 = DBg==CP.Channel_Concavity.Concavity_BasinIDs(i);
            bb = bwboundaries(t1);
            tx = [];
            ty = [];
            for j = 1:size(bb{1},1)
                tx = [tx;cutx(bb{1}(j,2))];
                ty = [ty;cuty(bb{1}(j,1))];
            end
            if sum(useBasinIDs == CP.Channel_Concavity.Concavity_BasinIDs(i)) > 0
                plot(tx./1000,ty./1000,'-r','linewidth',lw)
            else
                plot(tx./1000,ty./1000,'-k','linewidth',.5)
            end
        end
        caxis([0,2])
        cb = colorbar;
        ylabel(cb,'|Concavity|')
        title([tPre,'Plotted Basins'])
        savePlot(GI,gcf,'Concavity')
    end
catch
    warning('Could not plot Channel Concavity')
end
    
%% Chi Total Results
if sum(~isnan(CP.Chi.Total_Chi.Chi))>0
    try
        figure('Name','Total Channel Chi','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        % Channel Chi Map
        hold on
        imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false,'tickstokm',true)
        p0 = plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5);
        p1 = plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',1);
        [~,p2] = contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',1.5);
        plotc(kmTotChiS,CP.Chi.Total_Chi.Chi,'linewidth',5); colormap((crameri('lajolla')))

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'\chi (m)')
        setAxes(cb,0)
        box on
        title([tPre,sprintf('Total Channel \\chi\n(Best-Fitting M/N = %.2f)',CP.Chi.Total_Chi.MN)])
        legend([p1,p2,p0],{'Edifice Boundary','\chi Cutoff Elevation','Basin Boundaries'})

        savePlot(GI,gcf,'Total_Chi_Channels')
    catch
        warning('Could not plot Total Chi Channel Map')
    end

    try
        figure('Name','Projected Total Chi','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,CP.Chi.Total_Chi.Maximum_Chi,'colormap',(crameri('lajolla',255)),'tickstokm',true)
        p1 = plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',1);
        [~,p2] = contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',1.5);
        plot(kmTotChiS,'-k','linewidth',blw)

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Basin Max Total \chi (m)')
        setAxes(cb,0)
        box on
        title([tPre,sprintf('Basin Maximum Total \\chi\n(Best-Fitting M/N = %.2f)',CP.Chi.Total_Chi.MN)])
        legend([p1,p2],{'Edifice Boundary','\chi Cutoff Elevation'})

        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,CP.Chi.Total_Chi.Upstream_Chi,'colormap',(crameri('lajolla',255)),'tickstokm',true)
        p1 = plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k','linewidth',1);
        [~,p2] = contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',1.5);
        plot(kmTotChiS,'-k','linewidth',blw)

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        title([tPre,sprintf('Upstream-Projected Total \\chi\n(Best-Fitting M/N = %.2f)',CP.Chi.Total_Chi.MN)])
        setAxes(cb,0)
        box on
        legend([p1,p2],{'Edifice Boundary','\chi Cutoff Elevation'})
        savePlot(GI,gcf,'Total_Chi_Projected')
    catch
        warning('Could not plot Projected Total Chi Map')
    end
else
    warning('Cannot plot Projected Total Chi Map; no channels with chi')
end

%% Individua Chi Results
if isstruct(CP.Chi.Basin_Chi.Chi)
    basinMN = DP.Drainage_Basins.DB;
    basinMN.Z(:) = NaN;
    for i = 1:length(CP.Chi.Basin_Chi.Chi)
        basinMN.Z(DP.Drainage_Basins.DB.Z==CP.Chi.Basin_Chi.Chi(i).BasinID) = CP.Chi.Basin_Chi.Chi(i).MN;
    end

    try
        figure('Name','Individual Channel M/N and Chi','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        % M/N
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,basinMN,'colormap',flipud(crameri('lapaz',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'M/N')
        setAxes(cb,0)
        box on
        title([tPre,'Channel Best-Fitting M/N'])

        % Channel Chi Map
        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false,'tickstokm',true)
        plot(DP.Drainage_Basins.DBxy(:,1)./1000,DP.Drainage_Basins.DBxy(:,2)./1000,'-k','linewidth',.5)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        for i = 1:length(CP.Chi.Basin_Chi.Chi)
            tmpBasChiS = CP.Chi.Basin_Chi.Chi_S(i).S;
            tmpBasChiS.x = tmpBasChiS.x./1000;
            tmpBasChiS.y = tmpBasChiS.y./1000;
            plotc(tmpBasChiS,CP.Chi.Basin_Chi.Chi(i).Chi,'linewidth',5); 
        end
        contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',.5);
        colormap(gca,(crameri('lajolla')))

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'\chi (m)')
        setAxes(cb,0)
        box on
        title([tPre,'Individual Channel \chi'])
        savePlot(GI,gcf,'Individual_Chi_Channels')
    catch
        warning('Could not plot Individual Chi Channel Map')
    end

    try
        figure('Name','Projected Individual Chi','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,CP.Chi.Basin_Chi.Maximum_Chi,'colormap',(crameri('lajolla',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        for i = 1:length(CP.Chi.Basin_Chi.Chi_S)
            tmpBasChiS = CP.Chi.Basin_Chi.Chi_S(i).S;
            tmpBasChiS.x = tmpBasChiS.x./1000;
            tmpBasChiS.y = tmpBasChiS.y./1000;
            plot(tmpBasChiS,'-k','linewidth',blw)
        end
        contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',.5);

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Basin Max Individual \chi (m)')
        setAxes(cb,0)
        box on
        title([tPre,'Basin Maximum Individual \chi'])

        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,CP.Chi.Basin_Chi.Upstream_Chi,'colormap',(crameri('lajolla',255)),'tickstokm',true)
        plot(GP.boundaryXY(:,1)./1000,GP.boundaryXY(:,2)./1000,'-k')
        for i = 1:length(CP.Chi.Basin_Chi.Chi_S)
            tmpBasChiS = CP.Chi.Basin_Chi.Chi_S(i).S;
            tmpBasChiS.x = tmpBasChiS.x./1000;
            tmpBasChiS.y = tmpBasChiS.y./1000;
            plot(tmpBasChiS,'-k','linewidth',blw)
        end
        contour(cutX./1000,cutY./1000,cutZ,[1,1]*CP.Chi.Chi_Cutoff_Elevation,'-r','linewidth',.5);

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        title([tPre,'Upstream-Projected Individual \chi'])
        setAxes(cb,0)
        box on
        savePlot(GI,gcf,'Individual_Chi_Projected')
    catch
        warning('Could not plot Projected Individual Chi Map')
    end
else
    warning('Cannot plot Projected Indvidual Chi Map; no channels with chi')
end

%% Divide Ordering
if GI.inputs.Analyze_Divides
    if isempty(kmDVD)
        warning('Cannot plot Divide Ordring; no divides recognized')
    else
        try
            figure('Name','Topo Divide Ordering Metrics','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
            % Divide Asymmetry Map
            subplot(2,2,3)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            plotc(kmDVD,DV.Divide_AsymmetryIndex,'caxis',[0,1],'limit',[1000 inf]./1000)
            plot(kmS,'-k','linewidth',1)
            box on
    
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
    
            cb = colorbar;
            ylabel(cb,'Divide Asymmetry')
            setAxes(cb,0)
            colormap(gca,colormap((crameri('lajolla'))))
            title([tPre,'Divide Asymmetry Index'])
            
            % Divide Distance Map
            subplot(2,2,1)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            plotc(kmDVD,DV.Divide_Topo.DVD.distance./1e3,'limit',[1000 inf]./1000)
            plot(kmS,'-k','linewidth',1)
            box on
    
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
    
            cb = colorbar;
            ylabel(cb,'Divide Distance (km)')
            setAxes(cb,0)
            colormap(gca,flipud(crameri('davos')))
            title([tPre,'Divide Distance'])
    
            % Divide Elevation Map
            subplot(2,2,2)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            plotc(kmDVD,GRIDobj(cutX./1000,cutY./1000,cutZ),'limit',[1000 inf]./1000)
            plot(kmS,'-k','linewidth',1)
            box on
    
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
    
            cb = colorbar;
            ylabel(cb,'Divide Elevation (m)')
            setAxes(cb,0)
            colormap(gca,flipud(crameri('acton')));
            title([tPre,'Divide Elevation'])
    
            % Divide Order
            mo = nanmax(DV.Divide_Topo.DVD.order(:));
            mSpacing = (mo-1)/mo;
            mStart = 1+mSpacing/2;
            subplot(2,2,4)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            plotc(kmDVD,DV.Divide_Topo.DVD.order)
            plot(kmS,'-k','linewidth',1)
            box on
    
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
    
            cb = colorbar;
            ylabel(cb,'Divide Order')
            caxis([1,mo])
            if mo < 11
                set(cb,'ytick',mStart:mSpacing:mo)
                set(cb,'yticklabel',1:mo)
            else
                tmpMS = floor(mo/10);
                set(cb,'ytick',mStart:mSpacing*tmpMS:mo)
                
                set(cb,'yticklabel',1:tmpMS:mo);
            end
            setAxes(cb,0)
            colormap(gca,crameri('hawaii',mo))
            title([tPre,'Divide Order'])
    
            savePlot(GI,gcf,'Divide_Topo_Ordering')
        catch
            warning('Could not plot Divide Ordering')
        end
    end
end

%% Divide Asymmetry Statistics
if GI.inputs.Analyze_Divides
    if isempty(kmDVD)
        warning('Cannot plot Divide Asymmetry Statistics; no divides recognized')
    else
        try
            figure('Name','Divide Asymmetry Statistics','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
            colormap((crameri('lajolla')))
    
            % Divide Distance-Elevation-Asymmetry plot
            subplot(1,2,1)
            scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_AsymmetryIndex,'o','filled')
            box on
    
            xlabel('Divide Distance (km)')
            ylabel('Divide Elevation (m)')
            setAxes(gca,fs)
    
            cb = colorbar;
            ylabel(cb,'Asymmetry Index')
            setAxes(cb,0)
            caxis([0,1])
            title([tPre,'Divide Asymmetry Statistics'])
            axis square
            box on
    
            subplot(1,2,2)
            hold on
            p1 = plot(DV.DAI_Gamma.Bin_Midpoints,DV.DAI_Gamma.PDF,'o-k','linewidth',2,'markersize',14','markerfacecolor','k');
            p2 = plot(DV.DAI_Gamma_Corrected.Bin_Midpoints,DV.DAI_Gamma_Corrected.PDF,'o-r','linewidth',2,'markersize',12','markerfacecolor','r');
            xlim([0,1])
            setAxes(gca,fs);
            xlabel('DAI Bin')
            ylabel('PDF')
            legend([p1,p2],{sprintf('\\Gamma: %.2f',DV.DAI_Gamma.Gamma),sprintf('Corrected \\Gamma: %.2f',DV.DAI_Gamma_Corrected.Gamma)});
            title([tPre,'Divide Asymmetry Index PDF']);
            box on
            axis square
    
            savePlot(GI,gcf,'Divide_Asymmetry_Statistics')
        catch
            warning('Could not plot Divide Asymmetry Statistics')
        end
    end
end

%% Divide Chi Statistics
if GI.inputs.Analyze_Divides
    if isempty(kmDVD)
        warning('Cannot plot Divide Chi Difference Statistics; no divides recognized')
    else
        try
            figure('Name','Divide Chi Difference Statistics','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
            colormap((crameri('lajolla')))
    
            % Divide Total Chi-Difference Map
            subplot(2,2,1)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            try
                plotc(kmDVD,DV.Divide_TotalChiDifference)
            catch
                disp('Issues plotting Total Chi difference across divides')
            end
            plot(kmS,'-k','linewidth',1)
            box on
            
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
            
            cb = colorbar;
            ylabel(cb,'Divide Total \chi Difference')
            setAxes(cb,0)
            title([tPre,'Divide Total \chi Difference'])
    
            % Divide Individual Chi-Difference Map
            subplot(2,2,2)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            try
                plotc(kmDVD,DV.Divide_BasinChiDifference)
            catch
                disp('Issues plotting Individual Chi difference across divides')
            end
            plot(kmS,'-k','linewidth',1)
            box on
            
            xlabel('X (km)')
            ylabel('Y (km)')
            setAxes(gca,fs)
            
            cb = colorbar;
            ylabel(cb,'Divide Individual \chi Difference')
            setAxes(cb,0)
            title([tPre,'Divide Individual \chi Difference'])
    
            % Divide Distance-Elevation-Total Chi plot
            subplot(2,2,3)
            try
                scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_TotalChiDifference,'o','filled')
            catch
                disp('Issues plotting Total Chi difference across divides')
            end
            box on
            
            xlabel('Divide Distance (km)')
            ylabel('Divide Elevation (m)')
            setAxes(gca,fs)
            
            cb = colorbar;
            ylabel(cb,'Total \chi Difference (m)')
            setAxes(cb,0)
            title([tPre,'Total \chi Difference Divide Statistics II'])
            axis square
            box on
    
            % Divide Distance-Elevation-Individual Chi plot
            subplot(2,2,4)
            try
                scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_BasinChiDifference,'o','filled')
            catch
                disp('Issues plotting Individual Chi difference across divides')
            end
            box on
            
            xlabel('Divide Distance (km)')
            ylabel('Divide Elevation (m)')
            setAxes(gca,fs)
            
            cb = colorbar;
            ylabel(cb,'Individual \chi Difference (m)')
            setAxes(cb,0)
            title([tPre,'Individual \chi Difference Divide Statistics II'])
            axis square
            box on
            savePlot(GI,gcf,'Divide_ChiDifference_Statistics')
        catch
            warning('Could not plot Divide Chi Difference Statistics')
        end
    end
end
    
%% Junction Connectivity
if GI.inputs.Analyze_Divides
    if isempty(kmDVD)
        warning('Cannot plot Divide Junction Connectivity; no divides recognized')
    else
        try
            figure('Name','Junction Connectivity','units','normalized','outerposition',[0 0 1 1],'visible',visPlots)
            % Connectivity Map
            subplot(1,2,1)
            hold on
            imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false,'tickstokm',true)
            p1 = plot(kmDVD,'color',[0 0 0],'endpoints',false,'junction',false,'limit',[1000 Inf]./1000);
            p2 = scatter(DV.Junction_X_Y_C_Z_D_A(:,1)./1000,DV.Junction_X_Y_C_Z_D_A(:,2)./1000,50,DV.Junction_X_Y_C_Z_D_A(:,3),'filled','markeredgecolor','k');
            colormap(useCM)
            cb = colorbar;
            ylabel(cb,'Junction connectivity')
            title([tPre,'Junction Connectivity'])
            legend([p1(end),p2],{'Divide','Junction'})
            setAxes(gca,fs)
            cb = colorbar;
            ylabel(cb,'Junction Connectivity')
            setAxes(cb,0)
            colormap(gca,(crameri('bilbao')))
            xlabel('X (km)')
            ylabel('Y (km)')
    
            % Connectivity Plot
            subplot(1,2,2)
            plot(DV.Junction_X_Y_C_Z_D_A(:,3),DV.Junction_X_Y_C_Z_D_A(:,4),'ok','markerfacecolor','k')
            xlabel('Junction Connectivity')
            ylabel('Junction Elevation')
            title([tPre,'Junction Relationship'])
            axis tight
            axis square
            box on
            setAxes(gca,fs)
            savePlot(GI,gcf,'Junction_Connectivity')
        catch
            warning('Could not plot Divide Junction Connectivity')
        end
    end
end
end

function setAxes(h,fS)
    set(h,'xcolor','k')
    set(h,'ycolor','k')
    if fS ~= 0
        set(h,'fontsize',fS)
    end
end

function savePlot(GI,h,name)
    if ~isempty(GI.inputs.saveFigFolder)
        saveas(h,[GI.inputs.saveFigFolder,GI.inputs.figPrefix,name,'.fig']);
        exportgraphics(h, [GI.inputs.saveFigFolder,GI.inputs.figPrefix,name,'.png']);
        % saveas(h,[GI.inputs.saveFigFolder,GI.inputs.figPrefix,name,'.png']);
    end
end
