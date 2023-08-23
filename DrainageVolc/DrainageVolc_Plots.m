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
GP = res.GeneralParams;
TP = res.TopoParams;
DP = res.DrainageParams;
CP = res.ChannelParams;
DV = res.DivideParams;

useCM = viridis(255);

fs = 15;
lw = 3;
blw = 2;

%% Raw DEM
    figure('Name','Raw DEM','units','normalized','outerposition',[0 0 1 1])
    
    hold on
    imageschs(GP.DEM0,[],'colormap',demcmap(GP.DEM0.Z(:)))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k','linewidth',blw)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs);
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title('Raw Elevations')
    savePlot(GP,gcf,'Raw_Topography')
    
%% Topography and Hypsometry
try
    figure('Name','Elevation','units','normalized','outerposition',[0 0 1 1])
    % Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title('Clipped Topography')

    % Hypsometry
    subplot(1,2,2)
    plot(TP.ZHyps_Areas,TP.ZHyps_Vals,'-k','linewidth',lw)
    
    xlabel('Normalized Area')
    ylabel('Normalized Elevation')
    setAxes(gca,fs);
    axis tight
    axis square
    box on
    title('Topography Distribution (Hypsometry)')
    savePlot(GP,gcf,'Topography');
catch
    warning('Could not plot Topography and Hypsometry')
end
    
%% Slope and Distribution
try
    figure('Name','Slope','units','normalized','outerposition',[0 0 1 1])
    % Slope Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Slope,'colormap',useCM,'caxis',[0,min([40,max(TP.Slope.Z(:))])]);
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Slope (deg)')
    setAxes(gca,0)
    box on
    title('Slope')

    % Slope Distribution
    subplot(1,2,2)
    plot(TP.Slope_Hyps_Areas,TP.Slope_Hyps_Vals,'-k','linewidth',lw)
    
    xlabel('Normalized Area')
    ylabel('Slope (deg)')
    setAxes(gca,fs)
 
    axis tight
    axis square
    box on
    title('Slope Distribution')
    ylim([0,min([40,max(TP.Slope.Z(:))])])
    savePlot(GP,gcf,'Slope')
catch
    warning('Could not plot Slopes and Distribution')
end

%% Profile Curvature and Distribution
try
    i1 = find(TP.CProf_Hyps_Areas>=.01,1,'last');
    i2 = find(TP.CProf_Hyps_Areas<=.99,1,'first');

    rV = max(abs([TP.CProf_Hyps_Vals(i1),TP.CProf_Hyps_Vals(i2)]));

    figure('Name','Profile Curvature','units','normalized','outerposition',[0 0 1 1])
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Profile,'colormap',useCM,'caxis',[-rV,rV])
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Profile Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title('Profile Curvature')

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.CProf_Hyps_Areas,TP.CProf_Hyps_Vals,'-k','linewidth',lw)
    
    xlabel('Normalized Area')
    ylabel('Profile Curvature (m^{-1})')
    setAxes(gca,fs)
    axis tight
    axis square
    box on
    ylim([-rV,rV])
    title('Profile Curvature Distribution')
    savePlot(GP,gcf,'Profile_Curvature')
catch
    warning('Could not plot Profile Curvature and Distribution')
end

%% Planform Curvature and Distribution
try
    i1 = find(TP.CPlan_Hyps_Areas>=.01,1,'last');
    i2 = find(TP.CPlan_Hyps_Areas<=.99,1,'first');

    rV = max(abs([TP.CPlan_Hyps_Vals(i1),TP.CPlan_Hyps_Vals(i2)]));

    figure('Name','Planform Curvature','units','normalized','outerposition',[0 0 1 1])
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Planform,'colormap',useCM,'caxis',[-rV,rV])
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Planform Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title('Planform Curvature')

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.CPlan_Hyps_Areas,TP.CPlan_Hyps_Vals,'-k','linewidth',lw)
    
    xlabel('Normalized Area')
    ylabel('Planform Curvature (m^{-1})')
    setAxes(gca,fs)
    axis tight
    axis square
    box on
    ylim([-rV,rV])
    title('Planform Curvature Distribution')
    savePlot(GP,gcf,'Planform_Curvature')
catch
    warning('Could not plot Planform Curvature and Distribution')
end
    
%% Roughness and Distribution
try
    if length(GP.inputs.roughnessWindows)>4
        warning('Number of roughness windows exceeds figure max grids; plotting only four.')
    end
    
    maxSub = min([4,length(GP.inputs.roughnessWindows)]);
    
    figure('Name','Roughness Maps','units','normalized','outerposition',[0 0 1 1])
    % Roughness Maps
    for i = 1:maxSub
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.Roughness(i).Roughness_Grid,'colormap',useCM)
        plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Roughness (m)')
        setAxes(cb,0)
        title(sprintf('%.1f m Roughness',TP.Roughness(i).TrueWindowRes))
    end
    
    savePlot(GP,gcf,'Roughness_Maps')
    
    % Roughness Distributions
    figure('Name','Roughness Distributions','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    hold on
    legTitles = {};
    for i = 1:length(GP.inputs.roughnessWindows)
        plot(TP.Roughness(i).Hypsometry_Areas,TP.Roughness(i).Hypsometry_Values,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.Roughness(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Roughness')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'Location','northwest')
    box on
    title('Roughness Distribution')
    
    subplot(1,2,2)
    hold on
    legTitles = {};
    for i = 1:length(GP.inputs.roughnessWindows)
        plot(TP.Roughness(i).Hypsometry_Areas,TP.Roughness(i).Hypsometry_NormValues,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.Roughness(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Normalized Roughness')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'location','northwest')
    box on
    title('Normalized Roughness Distribution')
    savePlot(GP,gcf,'Roughness_Distributions')
catch
    warning('Could not plot Roughness and Distribution')
end

%% Windowed Slope Variance and Distribution
try
    if length(GP.inputs.slopeVarianceWindows)>4
        warning('Number of slope variance windows exceeds figure max grids; plotting only four.')
    end
    
    maxSub = min([4,length(GP.inputs.slopeVarianceWindows)]);
    
    figure('Name','Slope Variance Maps','units','normalized','outerposition',[0 0 1 1])
    % Slope Variance Maps
    for i = 1:maxSub
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.SlopeVariance_Windows(i).SlopeVariance_Grid,'colormap',useCM)
        plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Roughness (m)')
        setAxes(cb,0)
        title({sprintf('%.1f m',TP.SlopeVariance_Windows(i).TrueWindowRes);'Slope Variance'})
    end
    
    savePlot(GP,gcf,'SlopeVariance_Windowed_Maps')
    
    % Slope Variance Distributions
    figure('Name','Slope Variance Distributions','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    hold on
    legTitles = {};
    for i = 1:length(GP.inputs.slopeVarianceWindows)
        plot(TP.SlopeVariance_Windows(i).Hypsometry_Areas,TP.SlopeVariance_Windows(i).Hypsometry_Values,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.SlopeVariance_Windows(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Slope Variance')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'Location','northwest')
    box on
    title('Slope Variance Distribution')
    
    subplot(1,2,2)
    hold on
    legTitles = {};
    for i = 1:length(GP.inputs.slopeVarianceWindows)
        plot(TP.SlopeVariance_Windows(i).Hypsometry_Areas,TP.SlopeVariance_Windows(i).Hypsometry_NormValues,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.SlopeVariance_Windows(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Normalized Slope Variance')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'location','northwest')
    box on
    title('Normalized Slope Variance Distribution')
    savePlot(GP,gcf,'SlopeVariance_Windowed_Distributions')
catch
    warning('Could not plot Windowed Slope Variance and Distribution')
end
    
%% Elevation Slope Variance and Distribution
try
    figure('Name','Slope Variance Elevation & Basin','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    hold on
    yyaxis left
    plot(TP.SlopeVariance_Elevation.Values(:,5),TP.SlopeVariance_Elevation.Values(:,1),'-k','linewidth',2)
    setAxes(gca,fs)
    xlabel('Slope Variance')
    ylabel('Elevation (m)')
    axis tight
    axis square

    yyaxis right
    plot(TP.SlopeVariance_Elevation.Values(:,5),TP.SlopeVariance_Elevation.Values(:,2),'-k','linewidth',2)
    setAxes(gca,fs)
    xlabel('Slope Variance')
    ylabel('Normalized Elevation')
    axis tight
    axis square
    box on
    title({'Slope Variance by Contour';sprintf('Total Slope Variance = %.2f',TP.SlopeVariance_Total)})

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinSlopeVariance,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')

    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Slope Variance')
    setAxes(cb,0)
    title('Basin Slope Variance')
catch
    warning('Could not plot Elevation and Basin Slope Variance')
end
%% Shape Center Locations
try
    figure('Name','Shape Centers','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    p1 = plot(TP.TopographicCenter_XY(1),TP.TopographicCenter_XY(2),'ko','markerfacecolor',[1,1,1]*.8,'markersize',10);
    p2 = plot(TP.GeometricCenter_XY(1),TP.GeometricCenter_XY(2),'ks','markerfacecolor','r','markersize',10);
    p3 = plot(TP.VolumetricCenter_XY(1),TP.VolumetricCenter_XY(2),'k^','markerfacecolor','g','markersize',8);
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    legend([p1,p2,p3],{'Topographic Center','Geometric Center','Volumetric Center'})
    title('Volcano ''Center''')
    savePlot(GP,gcf,'Shape_Centers')
catch
    warning('Could not plot Shape Centers')
end
    
%% Drainage Area, Distance, and Basins
try
    figure('Name','Drainage Metrics','units','normalized','outerposition',[0 0 1 1])
    % Drainage Area Map
    subplot(2,2,1)
    hold on
    tmp = DP.A;
    tmp.Z = log10(tmp.Z*GP.inputs.dx^2);
    imageschs(GP.DEM,tmp,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,{'Log_{10} Cumulative';'Drainage Area (m^2)'})
    setAxes(cb,0)
    title('Cumulative Drainage Area')

    % Drainage Distance Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.D,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs);
    
    cb = colorbar;
    ylabel(cb,'Flow Distance (m)')
    setAxes(cb,0)
    title('Flow Distance')

    % Drainage Basins Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.DB,'colormap',jet(255))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Basin ID')
    setAxes(cb,0)
    title('Drainage Basins')

    % Drainage Basins Distribution
    subplot(2,2,4)
    hold on
    plot(DP.DB_Hyps_numDB,DP.DB_Hyps_Topo,'k','linewidth',2)
    plot(DP.DB_Hyps_numDB_Channelized,DP.DB_Hyps_Topo,'r','linewidth',2)
    box on
    set(gca,'xscale','log')
    setAxes(gca,fs)
    ylabel('Normalized Topography')
    xlabel('Number of Drainage Basins')
    axis tight
    axis square
    title('Basin Topographic Distribution')
    legend('All Basins','Channelized Basins','location','eastoutside')
    savePlot(GP,gcf,'Drainage_Metrics')
catch
    warning('Could not plot Drainage Values')
end
    
%% Hack's Law Plot - Basin Length
try
    tmpMap = DP.HackLawDeviation_BasinLength_Map;
    figure
    pcolor(tmpMap.Z);shading flat;
    bwrCM = colormap(bluewhitered(255));
    close

    HL = DP.HackLawDeviation_BasinLength;
    figure('Name','Hack''s Law Relationship (Basin Length)','units','normalized','outerposition',[0 0 1 1])
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
                tmpMap.Z(DP.DB.Z==HL(i,1)) = NaN;
            end
        end
    end
    tmp = plot(tmpRangeA,DP.HackLawFit_BasinLength(1)*tmpRangeA.^DP.HackLawFit_BasinLength(2),'-k','linewidth',2);
    pp = [pp;tmp];
    pL = [pL;{sprintf('Best Fit (L = %.2f A^{%.2f})',DP.HackLawFit_BasinLength(1),DP.HackLawFit_BasinLength(2))}];
    xlim(tmpRangeA)

    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,fs)
    
    xlabel('Basin Drainage Area (m^2)')
    ylabel('Basin Length (m)')
    axis tight
    axis square
    box on
    if GP.inputs.limitHacksLaw
        yy = ylim;
        tmp = plot([1,1]*CP.ChannelThreshold,yy,'--r','linewidth',2);
        pp = [pp;tmp];
        pL = [pL;{'Channelization Threshold'}];
    end
    legend(pp,pL,'location','northwest')
    title('Hack''s Law Relationship (Basin Length)')

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,tmpMap,'colormap',bwrCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title(sprintf('Basin Length Deviation from\nHack''s Law (Basin Length)'))
    savePlot(GP,gcf,'Hacks_Law_Relationship_BasinLength')
catch
    warning('Could not plot Basin Length Hack''s Law')
end

%% Hack's Law Plot - Basin Flow Length
try
    figure
    pcolor(DP.HackLawDeviation_FlowLength_Map.Z);shading flat;
    bwrCM = colormap(bluewhitered(255));
    close

    HL = DP.HackLawDeviation_FlowLength;
    figure('Name','Hack''s Law Relationship (Basin Flow Length)','units','normalized','outerposition',[0 0 1 1])
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
    tmp = plot(tmpRangeA,DP.HackLawFit_FlowLength(1)*tmpRangeA.^DP.HackLawFit_FlowLength(2),'-k','linewidth',2);
    pp = [pp;tmp];
    pL = [pL;{sprintf('Best Fit (L = %.2f A^{%.2f})',DP.HackLawFit_FlowLength(1),DP.HackLawFit_FlowLength(2))}];
    xlim(tmpRangeA)

    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,fs)

    xlabel('Basin Drainage Area (m^2)')
    ylabel('Max Basin Flow Length (m)')
    axis tight
    axis square
    box on
    if GP.inputs.limitHacksLaw
        yy = ylim;
        tmp = plot([1,1]*CP.ChannelThreshold,yy,'--r','linewidth',2);
        pp = [pp;tmp];
        pL = [pL;{'Channelization Threshold'}];
    end
    legend(pp,pL,'location','northwest')
    title('Hack''s Law Relationship (Basin Flow Length)')

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,DP.HackLawDeviation_FlowLength_Map,'colormap',bwrCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title(sprintf('Basin Length Deviation from\nHack''s Law (Flow Length)'))
    savePlot(GP,gcf,'Hacks_Law_Relationship_FlowLength')
catch
    warning('Could not plot Flow Length Hack''s Law')
end

%% Drainage Area - Slope Plots
try
    figure('Name','Slope-Area Plots','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    plot(DP.TopN_All_Area_Slope(:,1),tand(DP.TopN_All_Area_Slope(:,2)),'.k')
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    setAxes(gca,15)
    yy = ylim();
    xx = xlim();
    xlabel('Drainage Area (m^2)')
    ylabel('Slope')
    axis square
    box on
    if GP.inputs.basinTopN >= 1
        title(sprintf('All Slope-Area of\nLargest %d Basins',length(CP.Concavity_BasinIDs)));
    elseif GP.inputs.basinTopN > 0
        title(sprintf('All Slope-Area of\nLargest %d%% Basins\n(%d Basins)',GP.inputs.basinTopN*100,length(CP.Concavity_BasinIDs)));
    elseif GP.inputs.basinTopN < 0 && GP.inputs.basinTopN > -1
        title(sprintf('All Slope-Area of Basins\nin Upper %d%% of Topography\n(%d Basins)',abs(GP.inputs.basinTopN)*100,length(CP.Concavity_BasinIDs)));
    end

    subplot(1,2,2)
    hold on
    p1 = plot(DP.TopN_All_Area_Slope(:,1),tand(DP.TopN_All_Area_Slope(:,2)),'.','markeredgecolor',[.8 .8 .8]);
    p2 = plot(DP.TopN_Flow_Area_Slope(:,1),tand(DP.TopN_Flow_Area_Slope(:,2)),'.k');
    
    if ~isempty(DP.TopN_AreaThreshold_Rs)
        as = DP.TopN_Flow_Area_Slope;
        as(as(:,2)==0,:) = [];
        as(as(:,1)<DP.TopN_TransitionThreshold_A,:) = [];
        AS1 = as(as(:,1)<=DP.TopN_AreaThreshold_Rs(1),:);
        AS2 = as(as(:,1)>DP.TopN_AreaThreshold_Rs(1),:);
        
        mdl1 = fitlm(log10(AS1(:,1)),log10(tand(AS1(:,2))));
        mdl2 = fitlm(log10(AS2(:,1)),log10(tand(AS2(:,2))));
        
        aa1 = [DP.TopN_TransitionThreshold_A,DP.TopN_AreaThreshold_Rs(1)];
        aa2 = [DP.TopN_AreaThreshold_Rs(1),max(as(:,1))];
        p3 = plot(aa1,10^mdl1.Coefficients.Estimate(1)*aa1.^mdl1.Coefficients.Estimate(2),'-r','linewidth',3);
        p4 = plot(aa2,10^mdl2.Coefficients.Estimate(1)*aa2.^mdl2.Coefficients.Estimate(2),'-b','linewidth',3);
        
        ss1 = tand(as(:,2));
        
        p5 = plot([DP.TopN_AreaThreshold_Rs(1),DP.TopN_AreaThreshold_Rs(1)],[min(ss1),max(ss1)],'--k','linewidth',3);
        p6 = plot([1,1]*DP.TopN_TransitionThreshold_A,[min(ss1),max(ss1)],'--r','linewidth',3);
        
        ps = [p1,p2,p6,p3,p4,p5];
        pt = {'All Slope-Area','Flowpath Slope-Area','Transition Zone Start',...
        sprintf('Regression 1 (r^2 = %.2f)',DP.TopN_AreaThreshold_Rs(2)),...
        sprintf('Regression 2 (r^2 = %.2f)(M/N = %.2f)',DP.TopN_AreaThreshold_Rs(3),DP.TopN_MN),...
        sprintf('Area Threhold = %.2f km^2 (r'' = %.2f)',DP.TopN_AreaThreshold_Rs(1)./1e6,DP.TopN_AreaThreshold_Rs(4))};
        
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
%     title(sprintf('Largest %d Basin Flowpath Slope-Area',length(CP.Concavity_BasinIDs)));
    if GP.inputs.basinTopN >= 1
        title(sprintf('Flowpath Slope-Area of\nLargest %d Basins',length(CP.Concavity_BasinIDs)));
    elseif GP.inputs.basinTopN > 0
        title(sprintf('Flowpath Slope-Area of\nLargest %d%% Basins\n(%d Basins)',GP.inputs.basinTopN*100,length(CP.Concavity_BasinIDs)));
    elseif GP.inputs.basinTopN < 0 && GP.inputs.basinTopN > -1
        title(sprintf('Flowpath Slope-Area of Basins\nin Upper %d%% of Topography\n(%d Basins)',abs(GP.inputs.basinTopN)*100,length(CP.Concavity_BasinIDs)));
    end
    savePlot(GP,gcf,'Slope-Area')
catch
    warning('Could not plot Slope-Area')
end
    
%% Basin Statistics 1
try
    figure('Name','Basin Total Statistics 1','units','normalized','outerposition',[0 0 1 1])
    % Basin Height Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinHeights,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title('Basin Relief')
    
    % Basin Length Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinLengths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Length (m)')
    setAxes(cb,0)
    title('Basin Length')
    
    % Basin Width Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinWidths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title('Maximum Basin Width')

    % Basin Ellipticity Map
    tmp = DP.Basin_Statistics_Grids.BasinWidths;
    tmp.Z = tmp.Z./DP.Basin_Statistics_Grids.BasinLengths.Z;
    tmp.Z(tmp.Z>1) = NaN;
    subplot(2,2,4)
    hold on
    imageschs(GP.DEM,tmp,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Ellipticity')
    setAxes(cb,0)
    title('Basin Ellipticity')
        
    savePlot(GP,gcf,'Basin_Total_Statistics_1')
catch
    warning('Could not plot Basin Statistics 1')
end

%% Basin Statistics 2
try
    figure('Name','Basin Total Statistics 2','units','normalized','outerposition',[0 0 1 1])
    % Basin Hypsometry
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinHyps,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Hypsometry Integral')
    setAxes(cb,0)
    title('Basin Hypsometry')
    
    % Basin Mean Slope
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinSlopes,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    title('Basin Mean Slopes')
    
    % Drainage Area Map
    subplot(2,2,3)
    hold on
    tmp = DP.Basin_Statistics_Grids.DrainageArea;
    tmp.Z = log10(tmp.Z);
    imageschs(GP.DEM,tmp,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Log_{10} Drainage Area (m^2)')
    setAxes(cb,0)
    title('Basin Drainage Area')

    % Basin Orientation Map
    subplot(2,2,4)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinOrients,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Orientation (Azimuth)')
    setAxes(cb,0)
    title('Basin Orientation')

    savePlot(GP,gcf,'Basin_Total_Statistics_2')
catch
    warning('Could not plot Basin Statistics 2')
end
%% Basin Statistics 3
try
    figure('Name','Basin Total Statistics 3','units','normalized','outerposition',[0 0 1 1])
    % Basin Flow Length Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.FlowLengths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Flow Length (m)')
    setAxes(cb,0)
    title('Basin Flow Length')

    % Basin Euclidean Length Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinEuclideanLength,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Euclidean Length (m)')
    setAxes(cb,0)
    title('Basin Euclidean Length')
    cax = caxis;
    caxis([1,cax(2)])

    % Basin Sinuosity Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.FlowSinuosity,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Sinuosity')
    setAxes(cb,0)
    title('Basin Flow Sinuosity')
    cax = caxis;
    caxis([1,cax(2)])
    
    savePlot(GP,gcf,'Basin_Total_Statistics_3')
catch
    warning('Could not plot Basin Statistics 3')
end
%% Cross-Basin Statistics
try
    figure('Name','Basin Cross Statistics','units','normalized','outerposition',[0 0 1 1])
    % Basin Widths
    subplot(2,2,1)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    hold on
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    scatter(DP.Basin_Cross_Statistics(:,2),DP.Basin_Cross_Statistics(:,3),[],DP.Basin_Cross_Statistics(:,5),'o','filled')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title('Cross-Basin Widths')
    colormap(useCM)
    
    % Basin Relief
    subplot(2,2,2)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    hold on
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    scatter(DP.Basin_Cross_Statistics(:,2),DP.Basin_Cross_Statistics(:,3),[],DP.Basin_Cross_Statistics(:,6),'o','filled')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title('Cross-Basin Relief')

    % Basin Incision
    subplot(2,2,[3,4])
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    hold on
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    scatter(DP.Basin_Cross_Statistics(:,2),DP.Basin_Cross_Statistics(:,3),[],log10(DP.Basin_Cross_Statistics(:,7)),'o','filled')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    caxis([-3,0])
    cb = colorbar;
    ylabel(cb,'Log_{10} Incision Index')
    setAxes(cb,0)
    title({'Cross-Basin Incision Index'; '(Relief / Width)'})

    savePlot(GP,gcf,'Cross_Basin_Statistics')
catch
    warning('Could not plot Cross-Basin Statistics')
end

%% Drainage Basins Per Contour
try
    figure('Name','Basins Per Contour','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    hold on
    bCCpCLA = DP.Basin_Contour_ContourP_Count_Length_Area;
    bCCpCLA_C = DP.Basin_Contour_ContourP_Count_Length_Area_Channelized;
    p1 = plot(bCCpCLA(:,3),bCCpCLA(:,1),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw);
    p2 = plot(bCCpCLA_C(:,3),bCCpCLA_C(:,1),'-or','markerfacecolor','r','markersize',10,'linewidth',blw);
    xlabel('Number of Basins')
    ylabel('Contour (m)')
    setAxes(gca,fs)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])
    yyaxis right
    plot(bCCpCLA(:,3),bCCpCLA(:,2),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw)
    setAxes(gca,fs)
    if DP.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Basin_TopN+1),'--r','linewidth',lw)
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    legend([p1,p2],{'All Basins','Channelized Basins'})
    title('Number of Basins per Contour')

    subplot(1,2,2)
    hold on
    bCCpCLA = DP.Basin_Contour_ContourP_Count_Length_Area;
    p1 = plot(bCCpCLA(:,3)./bCCpCLA(:,4)*1000,bCCpCLA(:,1),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw);
    p2 = plot(bCCpCLA_C(:,3)./bCCpCLA_C(:,4)*1000,bCCpCLA_C(:,1),'-or','markerfacecolor','r','markersize',10,'linewidth',blw);
    xlabel('Basins per Contour Length (km^{-1})')
    setAxes(gca,fs)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])
    
    yyaxis right
    plot(bCCpCLA(:,3)./bCCpCLA(:,4)*1000,bCCpCLA(:,2),'-ok','markerfacecolor','k','markersize',10,'linewidth',blw)
    
    setAxes(gca,fs)
    title('Number of Basins per Contour Length')
    if DP.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Basin_TopN+1),'--r','linewidth',lw)
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    ylabel('Percent Relief')
    legend([p1,p2],{'All Basins','Channelized Basins'})
    savePlot(GP,gcf,'Basin_per_Contour')
catch
    warning('Could not plot Basin-Contour Relationships')
end

%% Radial Basin Count
try
    figure('Name','Radial Basin Count','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,DP.Radial_Analysis.Normalized_Radial_Distances,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(TP.TopographicCenter_XY(1),TP.TopographicCenter_XY(2),'ok','markerfacecolor','r')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Normalized Distance')
    setAxes(cb,0)
    title('Normalized Radial Distance')

    subplot(1,2,2)
    hold on
    plot(DP.Radial_Analysis.Basin_Count(:,2),DP.Radial_Analysis.Basin_Count(:,3),'-k','linewidth',lw)
    plot(DP.Radial_Analysis.Basin_Count(:,2),DP.Radial_Analysis.Basin_Count(:,4),'-r','linewidth',lw)
    
    xlabel('Normalized Radial Distance')
    ylabel('Basin Count')
    setAxes(gca,fs)
    legend('All Basins','Channelized Basins','location','northwest')
    title('Number of Radial Basins')
    box on
    axis square
    savePlot(GP,gcf,'Basin_per_RadialDistance')
catch
    warning('Could not plot Radial basin analysis')
end

%% Channels
try
    figure('Name','Channels','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k','linewidth',2)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title(sprintf('Drainage Channels (Threshold = %.2f km^2)',CP.ChannelThreshold./1e6))
    savePlot(GP,gcf,'Channel')
catch
    warning('Could not plot Channels')
end

%% Drainage Density Along Channels
try
    figure('Name','Channel Drainage Density','units','normalized','outerposition',[0 0 1 1])
    % Channel Drainage Density Map
    hold on
    imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false)
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plotc(CP.S,CP.DD,'linewidth',5); colormap(useCM)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Drainage Density (m^{-1})')
    setAxes(cb,0)
    box on
    title('Drainage Density Along Channels')

    savePlot(GP,gcf,'Drainage_Density_Channels')
catch
    warning('Could not plot Drainage Density Along Channels')
end

%% Basin Drainage Density
try
    figure('Name','Basin Drainage Density','units','normalized','outerposition',[0 0 1 1])
    % Basin Drainage Density Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,CP.BasinDD,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k','linewidth',blw)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Total Basin Drainage Density (m^{-1})')
    setAxes(cb,0)
    title(sprintf('Total Basin\nDrainage Density'))

    % Max Basin Drainage Density Map
    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,CP.MaxDD,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k','linewidth',blw)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Max Basin Drainage Density (m^{-1})')
    setAxes(cb,0)
    title(sprintf('Maximum Basin\nDrainage Density'))
    savePlot(GP,gcf,'Drainage_Density_Basins')
catch
    warning('Could not plot Basin Drainage Density')
end
    
%% Channel Concavity
try
    figure('Name','Concavity','units','normalized','outerposition',[0 0 1 1])
    
    if size(CP.Concavity_Streams,1) == 0
        useBasinIs = [];
    else
        nonNaNIs = 1:size(CP.Concavity_Stats);
        t = ones(1,size(CP.Concavity_Stats,1));
        for i = 1:size(CP.Concavity_Stats,1)
            t(i) = isnan(CP.Concavity_Stats(i).ks);
        end
        nonNaNIs(t==1) = [];

        if length(nonNaNIs) <= 4
            useBasinIs = length(nonNaNIs);
        else
            useBasinIs = randsample(nonNaNIs,4);
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
        cc = [CP.Concavity_Stats(useBasinIs(i)).a,CP.Concavity_Stats(useBasinIs(i)).g];
        cc(cc(:,2)<1e-3,:) = [];
        plot(cc(:,1),cc(:,2),'sk','markerfacecolor','r','markersize',10)
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xx = xlim();
        yy = ylim();
        plot(sortrows(cc(:,1)),CP.Concavity_Stats(i).ks(1)*sortrows(cc(:,1)).^CP.Concavity_Stats(i).theta,'--k','linewidth',lw);
        xlim(xx);
        ylim(yy);
        setAxes(gca,fs)
        title(sprintf('Basin %d\n\\theta = %.2f',CP.Concavity_BasinIDs(useBasinIs(i)),CP.Concavity_Stats(useBasinIs(i)).theta))
        box on
        if i > 2
            xlabel('Area (m^2)')
        end
        if i == 1 || i == 5
            ylabel('Slope')
        end
        axis square
        useBasinIDs = [useBasinIDs;CP.Concavity_BasinIDs(useBasinIs(i))];
    end
    
    % Concavity Basin Map
    subplot(2,4,[3,4,7,8])
    hold on
    try
        imageschs(GP.DEM,CP.Concavity_DEM,'colormap',useCM,'caxis',[0,2])
    catch
        warning('No Concavity values available - decrease the drainage area threshold.');
    end

    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k','linewidth',.5)
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    title('Drainage Basins')
    [DBg,~,~] = GRIDobj2mat(DP.DB);
    for i = 1:size(CP.Concavity_Streams,1)
        t1 = DBg==CP.Concavity_BasinIDs(i);
        bb = bwboundaries(t1);
        tx = [];
        ty = [];
        for j = 1:size(bb{1},1)
            tx = [tx;GP.cutX(bb{1}(j,2))];
            ty = [ty;GP.cutY(bb{1}(j,1))];
        end
        if sum(useBasinIDs == CP.Concavity_BasinIDs(i)) > 0
            plot(tx,ty,'-r','linewidth',lw)
        else
            plot(tx,ty,'-k','linewidth',.5)
        end
    end
    caxis([0,2])
    cb = colorbar;
    ylabel(cb,'|Concavity|')
    title('Plotted Basins')
    colormap(useCM)
    savePlot(GP,gcf,'Concavity')
catch
    warning('Could not plot Channel Concavity')
end
    
%% Chi Results
if ~isnan(CP.Chi)
    try
        figure('Name','Channel Chi','units','normalized','outerposition',[0 0 1 1])
        % Channel Chi Map
        hold on
        imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false)
        plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
        plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
        plotc(CP.chiS,CP.Chi,'linewidth',5); colormap(useCM)
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'\chi (m)')
        setAxes(cb,0)
        box on
        title(sprintf('Channel \\chi\n(Best-Fitting M/N = %.2f)',CP.BestFit_MN))
        savePlot(GP,gcf,'Chi_Channels')
    catch
        warning('Could not plot Chi Channel Map')
    end

    try
        figure('Name','Projected Chi','units','normalized','outerposition',[0 0 1 1])
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,CP.MaxChi,'colormap',useCM)
        plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
        plot(CP.chiS,'-k','linewidth',blw)
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Basin Max \chi (m)')
        setAxes(cb,0)
        box on
        title(sprintf('Basin Maximum \\chi\n(Best-Fitting M/N = %.2f)',CP.BestFit_MN))
    
        subplot(1,2,2)
        hold on
        imageschs(GP.DEM,CP.UpstreamChi,'colormap',useCM)
        plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
        plot(CP.chiS,'-k','linewidth',blw)
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        title(sprintf('Upstream-Projected \\chi\n(Best-Fitting M/N = %.2f)',CP.BestFit_MN))
        setAxes(cb,0)
        box on
        savePlot(GP,gcf,'Chi_Projected')
    catch
        warning('Could not plot Projected Chi Map')
    end
end

%% Divide Ordering
if GP.inputs.Analyze_Divides
    try
        figure('Name','Topo Divide Ordering Metrics','units','normalized','outerposition',[0 0 1 1])
        % Divide Distance Map
        subplot(2,2,1)
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
        plotc(DV.Divide_Topo.DVD,DV.Divide_Topo.DVD.distance./1e3,'limit',[1000 inf])
        plot(CP.S,'-k','linewidth',1)
        box on
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Divide Distance (km)')
        setAxes(cb,0)
        colormap(useCM)
        title('Divide Distance')
    
        % Divide Elevation Map
        subplot(2,2,2)
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
        plotc(DV.Divide_Topo.DVD,GP.DEM,'limit',[1000 inf])
        plot(CP.S,'-k','linewidth',1)
        box on
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Divide Elevation (m)')
        setAxes(cb,0)
        colormap(useCM)
        title('Divide Elevation')
    
        % Divide Asymmetry Map
        subplot(2,2,3)
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
        plotc(DV.Divide_Topo.DVD,DV.Divide_AsymmetryIndex,'caxis',[0,1],'limit',[1000 inf])
        plot(CP.S,'-k','linewidth',1)
        box on
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Divide Asymmetry')
        setAxes(cb,0)
        colormap(useCM)
        title('Divide Asymmetry Index')
    
        % Divide Chi-Difference Map
        subplot(2,2,4)
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
        try
            plotc(DV.Divide_Topo.DVD,DV.Divide_ChiDifference)
        catch
            disp('Issues plotting Chi difference across divides')
        end
        plot(CP.S,'-k','linewidth',1)
        box on
        
        xlabel('X (m)')
        ylabel('Y (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'Divide \chi Difference')
        setAxes(cb,0)
        colormap(useCM)
        title('Divide \chi Difference')
    
        savePlot(GP,gcf,'Divide_Topo_Ordering')
    catch
        warning('Could not plot Divide Ordering')
    end

    try
        figure('Name','Divide Statistics','units','normalized','outerposition',[0 0 1 1])
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
        title('Divide Topography Statistics I')
        savePlot(GP,gcf,'Divides_Topo_Ordering')
        axis square
        box on
    
        % Divide Distance-Elevation-Chi plot
        subplot(1,2,2)
        try
            scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_ChiDifference,'o','filled')
        catch
            disp('Issues plotting Chi difference across divides')
        end
        box on
        
        xlabel('Divide Distance (km)')
        ylabel('Divide Elevation (m)')
        setAxes(gca,fs)
        
        cb = colorbar;
        ylabel(cb,'\chi Difference (m)')
        setAxes(cb,0)
        title('Divide Topography Statistics II')
        colormap(useCM)
        axis square
        box on
        savePlot(GP,gcf,'Divide_Statistics')
    catch
        warning('Could not plot Divide Statistics')
    end
end
    
%% Junction Connectivity
if GP.inputs.Analyze_Divides
    try
        figure('Name','Junction Connectivity','units','normalized','outerposition',[0 0 1 1])
        % Connectivity Map
        subplot(1,2,1)
        hold on
        imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
        plot(DV.Divide_Topo.DVD,'color',[0 0 0],'endpoints',false,'junction',false)
        scatter(DV.Junction_X_Y_C_Z_D_A(:,1),DV.Junction_X_Y_C_Z_D_A(:,2),50,DV.Junction_X_Y_C_Z_D_A(:,3),'filled')
        colormap(useCM)
        cb = colorbar;
        ylabel(cb,'Junction connectivity')
        title('Junction Connectivity')
        legend('Divide','Junction')
        setAxes(gca,fs)
        cb = colorbar;
        ylabel(cb,'Junction Connectivity')
        setAxes(cb,0)
    
        % Connectivity Plot
        subplot(1,2,2)
        plot(DV.Junction_X_Y_C_Z_D_A(:,3),DV.Junction_X_Y_C_Z_D_A(:,4),'ok','markerfacecolor','k')
        xlabel('Junction Connectivity')
        ylabel('Junction Elevation')
        title('Junction Relationship')
        axis tight
        axis square
        box on
        setAxes(gca,fs)
        savePlot(GP,gcf,'Junction_Connectivity')
    catch
        warning('Could not plot Divide Junction Connectivity')
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

function savePlot(GP,h,name)
    if ~isempty(GP.inputs.saveFigFolder)
        saveas(h,[GP.inputs.saveFigFolder,GP.inputs.figPrefix,name,'.fig']);
        saveas(h,[GP.inputs.saveFigFolder,GP.inputs.figPrefix,name,'.png']);
    end
end