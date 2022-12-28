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

%% Raw DEM
    figure('Name','Raw DEM','units','normalized','outerposition',[0 0 1 1])
    
    hold on
    imageschs(GP.DEM0,[],'colormap',demcmap(GP.DEM0.Z(:)))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10);
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title('Raw Elevations')
    savePlot(GP,gcf,'Raw_Topography')
    
%% Topography and Hypsometry
    figure('Name','Elevation','units','normalized','outerposition',[0 0 1 1])
    % Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title('Clipped Topography')

    % Hypsometry
    subplot(1,2,2)
    plot(TP.ZHyps_Areas,TP.ZHyps_Vals,'-k','linewidth',2)
    
    xlabel('Normalized Area')
    ylabel('Normalized Elevation')
    setAxes(gca,10);
    axis tight
    axis square
    box on
    title('Topography Distribution')
    savePlot(GP,gcf,'Topography');
    
%% Slope and Distribution
    figure('Name','Slope','units','normalized','outerposition',[0 0 1 1])
    % Slope Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Slope,'colormap',useCM);
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Slope (deg)')
    setAxes(gca,0)
    box on
    title('Slope')

    % Slope Distribution
    subplot(1,2,2)
    plot(TP.Slope_Hyps_Areas,TP.Slope_Hyps_Vals,'-k','linewidth',2)
    
    xlabel('Normalized Area')
    ylabel('Slope (deg)')
    setAxes(gca,10)
 
    axis tight
    axis square
    box on
    title('Slope Distribution')
    savePlot(GP,gcf,'Slope')

%% Profile Curvature and Distribution
    rV = 1e-2;

    figure('Name','Profile Curvature','units','normalized','outerposition',[0 0 1 1])
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Profile,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Profile Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title('Profile Curvature')

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.CProf_Hyps_Areas,TP.CProf_Hyps_Vals,'-k','linewidth',2)
    
    xlabel('Normalized Area')
    ylabel('Profile Curvature (m^{-1})')
    setAxes(gca,10)
    axis tight
    axis square
    box on
    ylim([-.015,.015])
    title('Profile Curvature Distribution')
    savePlot(GP,gcf,'Profile_Curvature')

%% Planform Curvature and Distribution
    rV = 5e-2;

    figure('Name','Planform Curvature','units','normalized','outerposition',[0 0 1 1])
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Planform,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Planform Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title('Planform Curvature')

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.CPlan_Hyps_Areas,TP.CPlan_Hyps_Vals,'-k','linewidth',2)
    
    xlabel('Normalized Area')
    ylabel('Planform Curvature (m^{-1})')
    setAxes(gca,10)
    axis tight
    axis square
    box on
    ylim([-.1,.1])
    title('Planform Curvature Distribution')
    savePlot(GP,gcf,'Planform_Curvature')
    
%% Roughness and Distribution
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
        setAxes(gca,10)
        
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
        plot(TP.Roughness(i).Hypsometry_Areas,TP.Roughness(i).Hypsometry_Values,'-','linewidth',2)
        legTitles = [legTitles;sprintf('%.1f m',TP.Roughness(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Roughness')
    setAxes(gca,10)
    axis tight
    axis square
    legend(legTitles)
    box on
    title('Roughness Distribution')
    
    subplot(1,2,2)
    hold on
    legTitles = {};
    for i = 1:length(GP.inputs.roughnessWindows)
        plot(TP.Roughness(i).Hypsometry_Areas,TP.Roughness(i).Hypsometry_NormValues,'-','linewidth',2)
        legTitles = [legTitles;sprintf('%.1f m',TP.Roughness(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Normalized Roughness')
    setAxes(gca,10)
    axis tight
    axis square
    legend(legTitles)
    box on
    title('Normalized Roughness Distribution')
    savePlot(GP,gcf,'Roughness_Distributions')
    
%% Shape Center Locations
    figure('Name','Shape Centers','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    p1 = plot(TP.TopographicCenter_XY(1),TP.TopographicCenter_XY(2),'ko','markerfacecolor','k','markersize',10);
    p2 = plot(TP.GeometricCenter_XY(1),TP.GeometricCenter_XY(2),'ks','markerfacecolor','r','markersize',8);
    p3 = plot(TP.VolumetricCenter_XY(1),TP.VolumetricCenter_XY(2),'k^','markerfacecolor','g','markersize',6);
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    legend([p1,p2,p3],{'Topographic Center','Geometric Center','Volumetric Center'})
    title('Volcano ''Center''')
    savePlot(GP,gcf,'Shape_Centers')
    
%% Drainage Area, Distance, and Basins
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
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Log Cumulative Drainage Area (m^2)')
    setAxes(cb,0)
    title('Cumulative Drainage Area')

    % Drainage Distance Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,DP.D,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10);
    
    cb = colorbar;
    ylabel(cb,'Flow Distance (m)')
    setAxes(cb,0)
    title('Flow Distance')

    % Drainage Basins Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,DP.DB,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    title('Drainage Basins')

    % Drainage Basins Distribution
    subplot(2,2,4)
    plot(DP.DB_Hyps_numDB,DP.DB_Hyps_Topo,'k','linewidth',2)
    set(gca,'xscale','log')
    setAxes(gca,10)
    ylabel('Normalized Topography')
    xlabel('Number of Drainage Basins')
    axis tight
    axis square
    title('Basin Topographic Distribution')
    savePlot(GP,gcf,'Drainage_Metrics')
    
%% Hack's Law Plot - Basin Length
    figure('Name','Hack''s Law Relationship (Basin Length)','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)%
    hold on
    tmpRangeA = [min(DP.Basin_Statistics(:,2)),max(DP.Basin_Statistics(:,2))];
    plot(DP.Basin_Statistics(:,2),DP.Basin_Statistics(:,3),'rs','markerfacecolor','r')
    plot(tmpRangeA,DP.HackLawFit_BasinLength(1)*tmpRangeA.^DP.HackLawFit_BasinLength(2),'-k','linewidth',2)
    
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,10)
    
    xlabel('Basin Drainage Area (m^2)')
    ylabel('Basin Length (m)')
    legend('Volcano Basins',...
        sprintf('Best Fit (L = %.2f A^{%.2f})',DP.HackLawFit_BasinLength(1),DP.HackLawFit_BasinLength(2)),'location','northwest')
    axis tight
    axis square
    box on
    title('Hack''s Law Relationship (Basin Length)')

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,DP.HackLawDeviation_BasinLength,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title('Basin Length Deviation from Hack''s Law')
    savePlot(GP,gcf,'Hacks_Law_Relationship_BasinLength')

%% Hack's Law Plot - Basin Flow Length
    figure('Name','Hack''s Law Relationship (Basin Flow Length)','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)%
    hold on
    tmpRangeA = [min(DP.Basin_Statistics(:,2)),max(DP.Basin_Statistics(:,2))];
    plot(DP.Basin_Statistics(:,2),DP.Basin_Statistics(:,3),'rs','markerfacecolor','r')
    plot(tmpRangeA,DP.HackLawFit_FlowLength(1)*tmpRangeA.^DP.HackLawFit_FlowLength(2),'-k','linewidth',2)
    
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    setAxes(gca,10)
    
    xlabel('Basin Drainage Area (m^2)')
    ylabel('Max Basin Flow Length (m)')
    legend('Volcano Basins',...
        sprintf('Best Fit (L = %.2f A^{%.2f})',DP.HackLawFit_FlowLength(1),DP.HackLawFit_FlowLength(2)),'location','northwest')
    axis tight
    axis square
    box on
    title('Hack''s Law Relationship (Basin Flow Length)')

    subplot(1,2,2)
    hold on
    imageschs(GP.DEM,DP.HackLawDeviation_FlowLength,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Deviation from Power Law (m)')
    setAxes(gca,0)
    title('Basin Flow Length Deviation from Hack''s Law')
    savePlot(GP,gcf,'Hacks_Law_Relationship_FlowLength')
   


%% Drainage Area - Slope Plots
    figure('Name','Slope-Area Plots','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    plot(DP.TopN_All_Area_Slope(:,1),tand(DP.TopN_All_Area_Slope(:,2)),'.k')
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    setAxes(gca,10)
    yy = ylim();
    xx = xlim();
    xlabel('Drainage Area (m^2)')
    ylabel('Slope')
    axis square
    box on
    title(sprintf('Largest %d Basin All Slope-Area',length(CP.Concavity_BasinIDs)));
    
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
        p3 = plot(aa1,10^mdl1.Coefficients.Estimate(1)*aa1.^mdl1.Coefficients.Estimate(2),'-r','linewidth',2);
        p4 = plot(aa2,10^mdl2.Coefficients.Estimate(1)*aa2.^mdl2.Coefficients.Estimate(2),'-b','linewidth',2);
        
        ss1 = tand(as(:,2));
        
        p5 = plot([DP.TopN_AreaThreshold_Rs(1),DP.TopN_AreaThreshold_Rs(1)],[min(ss1),max(ss1)],'--k','linewidth',2);
        p6 = plot([1,1]*DP.TopN_TransitionThreshold_A,[min(ss1),max(ss1)],'--r','linewidth',2);
        
        ps = [p1,p2,p6,p3,p4,p5];
        pt = {'All A-S','Flowpath A-S','Transition Zone Start',...
        sprintf('Regression 1 (r^2 = %.2f)',DP.TopN_AreaThreshold_Rs(2)),...
        sprintf('Regression 2 (r^2 = %.2f)(M/N = %.2f)',DP.TopN_AreaThreshold_Rs(3),DP.TopN_MN),...
        sprintf('Area Threhold = %.2f km^2 (r'' = %.3f)',DP.TopN_AreaThreshold_Rs(1)./1e6,DP.TopN_AreaThreshold_Rs(4))};
        
    else
        ps = [p1,p2];
        pt = {'All A-S','Flowpath A-S'};
    end
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    setAxes(gca,10)
    xlim(xx)
    ylim(yy)
    xlabel('Drainage Area (m^2)')
    ylabel('Slope')
    axis square
    box on
    legend(ps,pt,'location','southwest')
    title(sprintf('Largest %d Basin Flowpath Slope-Area',length(CP.Concavity_BasinIDs)));
    savePlot(GP,gcf,'Slope-Area')
    
%% Basin Statistics 1
    figure('Name','Basin Total Statistics 1','units','normalized','outerposition',[0 0 1 1])
    % Basin Height Map
    subplot(2,3,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinHeights,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title('Basin Relief')
    
    % Basin Length Map
    subplot(2,3,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinLengths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Length (m)')
    setAxes(cb,0)
    title('Basin Length')
    
    % Basin Width Map
    subplot(2,3,3)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinWidths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title('Maximum Basin Width')

    % Basin Ellipticity Map
    tmp = DP.Basin_Statistics_Grids.BasinWidths;
    tmp.Z = tmp.Z./DP.Basin_Statistics_Grids.BasinLengths.Z;
    tmp.Z(tmp.Z>1) = NaN;
    subplot(2,3,4)
    hold on
    imageschs(GP.DEM,tmp,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Ellipticity')
    setAxes(cb,0)
    title('Basin Ellipticity')
    
    % Basin Orientation Map
    subplot(2,3,5)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinOrients,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Orientation (Azimuth)')
    setAxes(cb,0)
    title('Basin Orientation')
        
    savePlot(GP,gcf,'Basin_Total_Statistics_1')

%% Basin Statistics 2
    figure('Name','Basin Total Statistics 2','units','normalized','outerposition',[0 0 1 1])
    % Basin Hypsometry
    subplot(2,3,1)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinHyps,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Hypsometry Integral')
    setAxes(cb,0)
    title('Basin Hypsometry')
    
    % Basin Mean Slope
    subplot(2,3,2)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinSlopes,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    title('Basin Mean Slopes')
    
    % Basin Width Map
    subplot(2,3,3)
    hold on
    tmp = DP.Basin_Statistics_Grids.DrainageArea;
    tmp.Z = log10(tmp.Z);
    imageschs(GP.DEM,tmp,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Log_{10} Drainage Area (m^2)')
    setAxes(cb,0)
    title('Basin Drainage Area')

    % Basin Flow Length Map
    subplot(2,3,4)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.FlowLengths,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Flow Length (m)')
    setAxes(cb,0)
    title('Basin Flow Length')
    
    % Basin Sinuosity Map
    subplot(2,3,5)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.FlowSinuosity,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Sinuosity')
    setAxes(cb,0)
    title('Basin Flow Sinuosity')
    cax = caxis;
    caxis([1,cax(2)])

    % Basin Euclidean Length Map
    subplot(2,3,6)
    hold on
    imageschs(GP.DEM,DP.Basin_Statistics_Grids.BasinEuclideanLength,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Euclidean Length (m)')
    setAxes(cb,0)
    title('Basin Euclidean Length')
    cax = caxis;
    caxis([1,cax(2)])
    
    savePlot(GP,gcf,'Basin_Total_Statistics_2')
%% Cross-Basin Statistics
    figure('Name','Basin Cross Statistics','units','normalized','outerposition',[0 0 1 1])
    % Basin Widths
    subplot(1,2,1)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    hold on
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    scatter(DP.Basin_Cross_Statistics(:,2),DP.Basin_Cross_Statistics(:,3),[],DP.Basin_Cross_Statistics(:,5),'o','filled')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Width (m)')
    setAxes(cb,0)
    title('Cross-Basin Widths')
    colormap(useCM)
    
    % Basin Relief
    subplot(1,2,2)
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    hold on
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    scatter(DP.Basin_Cross_Statistics(:,2),DP.Basin_Cross_Statistics(:,3),[],DP.Basin_Cross_Statistics(:,6),'o','filled')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Relief (m)')
    setAxes(cb,0)
    title('Cross-Basin Relief')
    savePlot(GP,gcf,'Cross_Basin_Statistics')

%% Drainage Basins Per Contour
    figure('Name','Basins Per Contour','units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    bCCpCLA = DP.Basin_Contour_ContourP_Count_Length_Area;
    plot(bCCpCLA(:,3),bCCpCLA(:,1),'ok','markerfacecolor','k')
    xlabel('Number of Basins')
    ylabel('Contour (m)')
    setAxes(gca,10)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])
    yyaxis right
    plot(bCCpCLA(:,3),bCCpCLA(:,2),'ok','markerfacecolor','k')
    ylabel('Percent Relief')
    setAxes(gca,10)
    if DP.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Basin_TopN+1),'--r')
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    title('Number of Basins per Contour')

    subplot(1,2,2)
    bCCpCLA = DP.Basin_Contour_ContourP_Count_Length_Area;
    plot(bCCpCLA(:,3)./bCCpCLA(:,4),bCCpCLA(:,1),'ok','markerfacecolor','k')
    xlabel('Number of Basins per Contour Length')
    ylabel('Contour (m)')
    setAxes(gca,10)
    ylim([min(bCCpCLA(:,1)),max(bCCpCLA(:,1))])
    yyaxis right
    plot(bCCpCLA(:,3)./bCCpCLA(:,4),bCCpCLA(:,2),'ok','markerfacecolor','k')
    ylabel('Percent Relief')
    setAxes(gca,10)
    title('Number of Basins per Contour Length')
    if DP.Basin_TopN < 0
        xx = xlim;
        plot(xx,[1,1]*(DP.Basin_TopN+1),'--r')
    end
    ylim([min(bCCpCLA(:,2)),max(bCCpCLA(:,2))])
    savePlot(GP,gcf,'Basin_per_Contour')

%% Channels
    figure('Name','Channels','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z))
%     pcolor(GP.cutX,GP.cutY,GP.cutZ);shading flat; axis image; colormap(demcmap(GP.cutZ));
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k','linewidth',2)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title(sprintf('Drainage Channels (Threshold = %.2f km^2)',CP.ChannelThreshold./1e6))
    savePlot(GP,gcf,'Channel')

%% Drainage Density
    figure('Name','Drainage Density','units','normalized','outerposition',[0 0 1 1])
    % Channel Drainage Density Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false)
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plotc(CP.S,CP.DD,'linewidth',2); colormap(useCM)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Drainage Density (m^{-1})')
    setAxes(cb,0)
    box on
    title('Drainage Density')

    % Basin Drainage Density Map
    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,CP.BasinDD,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Basin Drainage Density (m^{-1})')
    setAxes(cb,0)
    title('Basin Drainage Density')

    % Max Basin Drainage Density Map
    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,CP.MaxDD,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.S,'-k')
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Max Drainage Density in Basin (m^{-1})')
    setAxes(cb,0)
    title('Maximum Drainage Density in Basin')
    savePlot(GP,gcf,'Drainage_Density')
    
%% Channel Concavity
    figure('Name','Concavity','units','normalized','outerposition',[0 0 1 1])
    if size(CP.Concavity_Streams,1) <= 8
        useBasinIs = 1:size(CP.Concavity_Streams,1);
    else
        useBasinIs = randsample(size(CP.Concavity_Streams,1),8);
    end
    useBasinIDs = [];
    % Concavity Plot
    for i = 1:min([size(CP.Concavity_Streams,1),8])
        if i > 4
            useI = i+2;
        else
            useI = i;
        end
        subplot(2,6,useI)
        hold on
        cc = [CP.Concavity_Stats(useBasinIs(i)).a,CP.Concavity_Stats(useBasinIs(i)).g];
        cc(cc(:,2)<1e-3,:) = [];
        plot(cc(:,1),cc(:,2),'sk','markerfacecolor','r')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xx = xlim();
        yy = ylim();
%         plot(xx,CP.Concavity_Stats(i).ks(1)*xx.^CP.Concavity_Stats(i).theta,'--k','linewidth',1.5);
        plot(sortrows(cc(:,1)),CP.Concavity_Stats(i).ks(1)*sortrows(cc(:,1)).^CP.Concavity_Stats(i).theta,'--k','linewidth',1.5);
        xlim(xx);
        ylim(yy);
        setAxes(gca,10)
        title(sprintf('Basin %d, \\theta = %.2f',CP.Concavity_BasinIDs(useBasinIs(i)),CP.Concavity_Stats(useBasinIs(i)).theta))
        box on
        xlabel('Drainage Area (m^2)')
        ylabel('Slope')
        axis square
        useBasinIDs = [useBasinIDs;CP.Concavity_BasinIDs(useBasinIs(i))];
    end
    
    % Concavity Basin Map
    subplot(2,6,[5,6,11,12])
    hold on
    try
        imageschs(GP.DEM,CP.Concavity_DEM,'colormap',useCM)
    catch
        warning('No Concavity values available - decrease the drainage area threshold.');
    end

    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k','linewidth',.5)
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
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
            plot(tx,ty,'-r','linewidth',2)
        else
            plot(tx,ty,'-k','linewidth',.5)
        end
    end
    caxis([0,1])
    cb = colorbar;
    ylabel(cb,'|Concavity|')
    title('Plotted Basins')
    colormap(useCM)
    savePlot(GP,gcf,'Concavity')
    
%% Chi Results
if ~isnan(CP.Chi)
    figure('Name','Chi','units','normalized','outerposition',[0 0 1 1])
    % Channel Chi Map
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',[1 1 1],'colorbar',false)
    plot(DP.DBxy(:,1),DP.DBxy(:,2),'-k','linewidth',.5)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plotc(CP.chiS,CP.Chi,'linewidth',2); colormap(useCM)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'\chi (m)')
    setAxes(cb,0)
    box on
    title(sprintf('Channel \\chi (Best-Fitting M/N = %.2f)',CP.BestFit_MN))

    subplot(2,2,2)
    hold on
    imageschs(GP.DEM,CP.MaxChi,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.chiS,'-k','linewidth',1)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Basin Max \chi (m)')
    setAxes(cb,0)
    box on
    title('Basin Maximum \chi')
    savePlot(GP,gcf,'Chi')

    subplot(2,2,3)
    hold on
    imageschs(GP.DEM,CP.UpstreamChi,'colormap',useCM)
    plot(GP.boundaryXY(:,1),GP.boundaryXY(:,2),'-k')
    plot(CP.chiS,'-k','linewidth',1)
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Upstream-Projected \chi values (m)')
    setAxes(cb,0)
    box on
    title('Upstream-Project \chi')
    savePlot(GP,gcf,'Chi')
end

%% Divide Ordering
if GP.inputs.Analyze_Divides
    figure('Name','Topo Divide Ordering Metrics','units','normalized','outerposition',[0 0 1 1])
    % Divide Distance Map
    subplot(2,3,1)
    hold on
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    plotc(DV.Divide_Topo.DVD,DV.Divide_Topo.DVD.distance./1e3,'limit',[1000 inf])
    plot(CP.S,'-k','linewidth',1)
%     plot(DV.Divide_Topo.X,DV.Divide_Topo.Y,'*k')
    box on
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Divide Distance (km)')
    setAxes(cb,0)
    colormap(useCM)
    title('Divide Distance')

    % Divide Elevation Map
    subplot(2,3,2)
    hold on
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    plotc(DV.Divide_Topo.DVD,GP.DEM,'limit',[1000 inf])
    plot(CP.S,'-k','linewidth',1)
    box on
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Divide Elevation (m)')
    setAxes(cb,0)
    colormap(useCM)
    title('Divide Elevation')

    % Divide Asymmetry Map
    subplot(2,3,3)
    hold on
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    plotc(DV.Divide_Topo.DVD,DV.Divide_AsymmetryIndex,'caxis',[0,1],'limit',[1000 inf])
    plot(CP.S,'-k','linewidth',1)
    box on
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Divide Asymmetry')
    setAxes(cb,0)
    colormap(useCM)
    title('Divide Asymmetry Index')

    % Divide Chi-Difference Map
    subplot(2,3,4)
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
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Divide \chi Difference')
    setAxes(cb,0)
    colormap(useCM)
    title('Divide \chi Difference')

    % Divide Distance-Elevation-Asymmetry plot
    subplot(2,3,5)
    scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_AsymmetryIndex,'o','filled')
    box on
    
    xlabel('Divide Distance (km)')
    ylabel('Divide Elevation (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'Asymmetry Index')
    setAxes(cb,0)
    caxis([0,1])
    title('Divide Topography Statistics I')
    savePlot(GP,gcf,'Divides_Topo_Ordering')

    % Divide Distance-Elevation-Chi plot
    subplot(2,3,6)
    try
        scatter(DV.Divide_Topo.DVD.distance./1000,getvalue(DV.Divide_Topo.DVD,DV.VerticalDistance,'min'),[],DV.Divide_ChiDifference,'o','filled')
    catch
        disp('Issues plotting Chi difference across divides')
    end
    box on
    
    xlabel('Divide Distance (km)')
    ylabel('Divide Elevation (m)')
    setAxes(gca,10)
    
    cb = colorbar;
    ylabel(cb,'\chi Difference')
    setAxes(cb,0)
    title('Divide Topography Statistics II')
    savePlot(GP,gcf,'Divides_Topo_Ordering')
end
    
%% Junction Connectivity
if GP.inputs.Analyze_Divides
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
    setAxes(gca,10)

    % Connectivity Plot
    subplot(1,2,2)
    plot(DV.Junction_X_Y_C_Z_D_A(:,3),DV.Junction_X_Y_C_Z_D_A(:,4),'ok','markerfacecolor','k')
    xlabel('Junction Connectivity')
    ylabel('Junction Elevation')
    title('Junction Relationship')
    axis tight
    axis square
    box on
    setAxes(gca,10)
    savePlot(GP,gcf,'Junction_Connectivity')
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