function MorVolc_Plots(res)
%%
% Name: MorVolc_Plots
% Date: 06/30/2021 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to plot results of the MorVolc analysis.
%
% Input:
%   res: Result package from MorVolc.

%% Unpack Structure
GI = res.GeneralParams;
GP = res.GeographicParams;
OP = res.OrientParams;
TP = res.TopoParams;
SP = res.SlopeParams;
ShP = res.ShapeParams;
SiP = res.SizeParams;
PP = res.PeakParams;
RP = res.RoughnessParams;
WsvP = res.WindowedSlopeVarianceParams;

useCM = viridis(255);

[Z,x,y] = GRIDobj2mat(GP.DEM);
[S,~,~] = GRIDobj2mat(GP.Slope);
[X,Y] = meshgrid(x,y);

figPre = GI.inputs.figTitlePrefix;
if ~isempty(figPre) && ~strcmp(figPre(end),' ')
    figPre = [figPre,' '];
end

fs = 15;
lw = 3;
blw = 2;

%% Hillshade DEM
try
    figure('Name','DEM Hillshade','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    imageschs(GP.DEM,[],'colormap',[1 1 1]*.9,'colorbar',false,'tickstokm',1)
    p1 = plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k','linewidth',blw);
    p2 = plot(GP.SummitXYZ(:,1)./1000,GP.SummitXYZ(:,2)./1000,'-r','linewidth',blw);
    p4 = plot(GP.Lower_Flank_XY(:,1)./1000,GP.Lower_Flank_XY(:,2)./1000,'-b','linewidth',blw);
    if ~isempty(GP.maskXY)
        for i = 1:length(GP.maskXY)
            p3 = fill(GP.maskXY{i}(:,1)./1000,GP.maskXY{i}(:,2)./1000,'w');
        end
        allP = [p1,p2,p4,p3];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary','Mask Region'};
    else
        allP = [p1,p2,p4];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary'};
    end
    if ~isempty(GP.craterXYZ)
        for i = 1:length(GP.craterXYZ)
            pp = fill(GP.craterXYZ{i}(:,1)./1000,GP.craterXYZ{i}(:,2)./1000,'r','facealpha',.5);
        end
        allP = [allP,pp];
        allPT = [allPT,{'Crater Region'}];
    end

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)
    title([figPre,'Edifice Hillshade'])
    legend(allP,allPT);
    savePlot(GI,gcf,'Hillshade_DEM')
catch
    warning('Could not plot Hillshade DEM')
end

%% 3D DEM
try
    figure('Name','3D DEM','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    surf(X./1000,Y./1000,Z); shading flat; colormap(demcmap(GP.DEM.Z)); camlight
    hold on
    p1 = plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,GP.boundaryXYZ(:,3),'-k','linewidth',blw);
    p2 = plot3(GP.SummitXYZ(:,1)./1000,GP.SummitXYZ(:,2)./1000,GP.SummitXYZ(:,3),'-r','linewidth',blw);
    p4 = plot3(GP.Lower_Flank_XY(:,1)./1000,GP.Lower_Flank_XY(:,2)./1000,GP.Lower_Flank_XY(:,3),'-b','linewidth',blw);
    title([figPre,'3D Topography'])
    legend([p1,p2,p4],{'Edifice Boundary','Summit Boundary','Lower Flank Boundary'});
    box on

    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    grid on
    axis tight
    savePlot(GI,gcf,'3D_DEM')
catch
    warning('Could not plot 3D DEM')
end

%% Topography and Hypsometry
try
    figure('Name','Elevation','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1)
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title([figPre,'Clipped Topography'])

    % Hypsometry
    subplot(1,2,2)
    plot(TP.Hypsometry.Elevation_Hyps_Areas_Values(:,1),TP.Hypsometry.Elevation_Hyps_Areas_Values(:,2),'-k','linewidth',lw)

    xlabel('Normalized Area')
    ylabel('Normalized Elevation')
    setAxes(gca,fs);
    axis tight
    axis square
    box on
    title([figPre,'Topography Distribution (Hypsometry)'])
    savePlot(GI,gcf,'Topography');
catch
    warning('Could not plot Topography and Hypsometry')
end

%% Slope and Distribution
try
    figure('Name','Slope','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Slope Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,GP.Slope,'colormap',flipud(crameri('lajolla',255)),'caxis',[0,min([40,max(GP.Slope.Z(:))])],'tickstokm',1);
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Slope (deg)')
    setAxes(gca,0)
    box on
    title([figPre,'Slope'])

    % Slope Distribution
    subplot(1,2,2)
    plot(TP.Hypsometry.Slope_Hyps_Areas_Values(:,1),TP.Hypsometry.Slope_Hyps_Areas_Values(:,2),'-k','linewidth',lw)

    xlabel('Normalized Area')
    ylabel('Slope (deg)')
    setAxes(gca,fs)

    axis tight
    axis square
    box on
    title([figPre,'Slope Distribution'])
    ylim([0,min([40,max(GP.Slope.Z(:))])])
    savePlot(GI,gcf,'Slope')
catch
    warning('Could not plot Slopes and Distribution')
end

%% 3D Slope
try
    figure('Name','3D Slope','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    surf(X./1000,Y./1000,Z,S); shading flat; colormap(flipud(crameri('lajolla',255))); camlight
    hold on
    p1 = plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,GP.boundaryXYZ(:,3),'-k','linewidth',2);
    p2 = plot3(GP.SummitXYZ(:,1)./1000,GP.SummitXYZ(:,2)./1000,GP.SummitXYZ(:,3),'-r','linewidth',2);
    p4 = plot3(GP.Lower_Flank_XY(:,1)./1000,GP.Lower_Flank_XY(:,2)./1000,GP.Lower_Flank_XY(:,3),'--k','linewidth',2);
    title([figPre,'3D Slope'])
    legend([p1,p2,p4],{'Edifice Boundary','Summit Boundary','Lower Flank Boundary'});

    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    caxis([0,40])
    box on
    grid on
    axis tight
    savePlot(GI,gcf,'3D_Slope')
catch
    warning('Could not plot 3D Slopes')
end

%% Profile Curvature and Distribution
try
    i1 = find(TP.Hypsometry.ProfileC_Hyps_Areas_Values(:,1)>=.01,1,'last');
    i2 = find(TP.Hypsometry.ProfileC_Hyps_Areas_Values(:,1)<=.99,1,'first');

    rV = max(abs([TP.Hypsometry.ProfileC_Hyps_Areas_Values(i1,2),TP.Hypsometry.ProfileC_Hyps_Areas_Values(i2,2)]));

    figure('Name','Profile Curvature','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Grids.Profile,'colormap',crameri('vik',255),'caxis',[-rV,rV],'tickstokm',1)
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Profile Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title([figPre,'Profile Curvature'])

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.Hypsometry.ProfileC_Hyps_Areas_Values(:,1),TP.Hypsometry.ProfileC_Hyps_Areas_Values(:,2),'-k','linewidth',lw)

    xlabel('Normalized Area')
    ylabel('Profile Curvature (m^{-1})')
    setAxes(gca,fs)
    axis tight
    axis square
    box on
    ylim([-rV,rV])
    title([figPre,'Profile Curvature Distribution'])
    savePlot(GI,gcf,'Profile_Curvature')
catch
    warning('Could not plot Profile Curvature and Distribution')
end

%% Planform Curvature and Distribution
try
    i1 = find(TP.Hypsometry.PlanformC_Hyps_Areas_Values(:,1)>=.01,1,'last');
    i2 = find(TP.Hypsometry.PlanformC_Hyps_Areas_Values(:,1)<=.99,1,'first');

    rV = max(abs([TP.Hypsometry.PlanformC_Hyps_Areas_Values(i1,2),TP.Hypsometry.PlanformC_Hyps_Areas_Values(i2,2)]));

    figure('Name','Planform Curvature','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Curvature Map
    subplot(1,2,1)
    hold on
    imageschs(GP.DEM,TP.Curvature_Grids.Planform,'colormap',crameri('bam',255),'caxis',[-rV,rV],'tickstokm',1)
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Planform Curvature (m^{-1})')
    setAxes(cb,0)
    caxis([-rV,rV])
    title([figPre,'Planform Curvature'])

    % Curvature Distribution
    subplot(1,2,2)
    hold on
    plot(TP.Hypsometry.PlanformC_Hyps_Areas_Values(:,1),TP.Hypsometry.PlanformC_Hyps_Areas_Values(:,2),'-k','linewidth',lw)

    xlabel('Normalized Area')
    ylabel('Planform Curvature (m^{-1})')
    setAxes(gca,fs)
    axis tight
    axis square
    box on
    ylim([-rV,rV])
    title([figPre,'Planform Curvature Distribution'])
    savePlot(GI,gcf,'Planform_Curvature')
catch
    warning('Could not plot Planform Curvature and Distribution')
end

%% Roughness and Distribution
try
    if length(GI.inputs.roughnessWindows)>4
        warning('Number of roughness windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.roughnessWindows)]);

    figure('Name','Roughness Maps','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Roughness Maps
    for i = 1:maxSub
        gM = max(abs(TP.Roughness(i).Roughness_Grid.Z(:)));
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.Roughness(i).Roughness_Grid,'colormap',flipud(crameri('roma',255)),'caxis',[-1,1]*gM,'tickstokm',1)
        plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Roughness (m)')
        setAxes(cb,0)
        title([figPre,sprintf('%.1f m Roughness',TP.Roughness(i).TrueWindowRes)])
    end

    savePlot(GI,gcf,'Roughness_Maps')

    % Roughness Distributions
    figure('Name','Roughness Distributions','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    subplot(1,2,1)
    hold on
    legTitles = {};
    for i = 1:length(GI.inputs.roughnessWindows)
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
    title([figPre,'Roughness Distribution'])

    subplot(1,2,2)
    hold on
    legTitles = {};
    for i = 1:length(GI.inputs.roughnessWindows)
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
    title([figPre,'Normalized Roughness Distribution'])
    savePlot(GI,gcf,'Roughness_Distributions')
catch
    warning('Could not plot Roughness and Distribution')
end

%% Windowed Slope Variance and Distribution
try
    if length(GI.inputs.slopeVarianceWindows)>4
        warning('Number of slope variance windows exceeds figure max grids; plotting only four.')
    end

    maxSub = min([4,length(GI.inputs.slopeVarianceWindows)]);

    figure('Name','Slope Variance Maps','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    % Slope Variance Maps
    for i = 1:maxSub
        gM = max(abs(TP.SlopeVariance.Windows(i).SlopeVariance_Grid.Z(:)));
        subplot(2,2,i)
        hold on
        imageschs(GP.DEM,TP.SlopeVariance.Windows(i).SlopeVariance_Grid,'colormap',flipud(crameri('roma',255)),'caxis',[-1,1]*gM,'tickstokm',1)
        plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')

        xlabel('X (km)')
        ylabel('Y (km)')
        setAxes(gca,fs)

        cb = colorbar;
        ylabel(cb,'Slope Variance')
        setAxes(cb,0)
        title({[figPre,sprintf('%.1f m',TP.SlopeVariance.Windows(i).TrueWindowRes)];'Slope Variance'})
    end

    savePlot(GI,gcf,'SlopeVariance_Windowed_Maps')

    % Slope Variance Distributions
    figure('Name','Slope Variance Distributions','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    subplot(1,2,1)
    hold on
    legTitles = {};
    for i = 1:length(GI.inputs.slopeVarianceWindows)
        plot(TP.SlopeVariance.Windows(i).Hypsometry_Areas,TP.SlopeVariance.Windows(i).Hypsometry_Values,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.SlopeVariance.Windows(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Slope Variance')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'Location','northwest')
    box on
    title([figPre,'Slope Variance Distribution'])

    subplot(1,2,2)
    hold on
    legTitles = {};
    for i = 1:length(GI.inputs.slopeVarianceWindows)
        plot(TP.SlopeVariance.Windows(i).Hypsometry_Areas,TP.SlopeVariance.Windows(i).Hypsometry_NormValues,'-','linewidth',lw)
        legTitles = [legTitles;sprintf('%.1f m',TP.SlopeVariance.Windows(i).TrueWindowRes)];
    end
    xlabel('Normalized Area')
    ylabel('Normalized Slope Variance')
    setAxes(gca,fs)
    axis tight
    axis square
    legend(legTitles,'location','northwest')
    box on
    title([figPre,'Normalized Slope Variance Distribution'])
    savePlot(GI,gcf,'SlopeVariance_Windowed_Distributions')
catch
    warning('Could not plot Windowed Slope Variance and Distribution')
end

%% Shape Center Locations
try
    figure('Name','Shape Centers','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1)
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k')
    p1 = plot(TP.Centers.Topographic_XY(1)./1000,TP.Centers.Topographic_XY(2)./1000,'ko','markerfacecolor',[1,1,1]*.8,'markersize',12);
    p2 = plot(TP.Centers.Geometric_XY(1)./1000,TP.Centers.Geometric_XY(2)./1000,'ks','markerfacecolor','r','markersize',12);
    p3 = plot(TP.Centers.Volumetric_XY(1)./1000,TP.Centers.Volumetric_XY(2)./1000,'k^','markerfacecolor','g','markersize',10);

    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    legend([p1,p2,p3],{'Topographic Center','Geometric Center','Volumetric Center'})
    title([figPre,'Volcano ''Center'''])
    savePlot(GI,gcf,'Shape_Centers')
catch
    warning('Could not plot Shape Centers')
end

%% Basal Surface Points 
try
    scarseInd = 5;
    tX = X(1:scarseInd:end,1:scarseInd:end);
    tY = Y(1:scarseInd:end,1:scarseInd:end);
    tZ = Z(1:scarseInd:end,1:scarseInd:end);
    pointsXYZ = [tX(:),tY(:),tZ(:)];
    pointsXYZ(isnan(pointsXYZ(:,3)),:) = [];

    edRel = range(pointsXYZ(:,3));
    edXYZ = pointsXYZ;
    edXYZ(edXYZ(:,3)>min(edXYZ(:,3))+edRel*.4,:) = [];

    figure('Name','Basal Surface','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    scarseInd = scarseInd*2;
    legendHit = 0;
    if GI.inputs.interpSurfaces.Natural
        [tmpIZ,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Natural);
        subplot(2,2,1)
        hold on
        plot3(edXYZ(:,1)./1000,edXYZ(:,2)./1000,edXYZ(:,3),'.k')
        plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,GP.boundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(X(1:scarseInd:end,1:scarseInd:end)./1000,Y(1:scarseInd:end,1:scarseInd:end)./1000,tmpIZ(1:scarseInd:end,1:scarseInd:end),'.r')

        xlabel('X (km)')
        ylabel('Y (km)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
        legendHit = 1;
        title([figPre,'Natural Interpolation'])
        view(20,20);
    end

    if GI.inputs.interpSurfaces.IDW
        [tmpIZ,~,~] = GRIDobj2mat(GP.Basal_Surfaces.IDW);
        subplot(2,2,2)
        hold on
        plot3(edXYZ(:,1)./1000,edXYZ(:,2)./1000,edXYZ(:,3),'.k')
        plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,GP.boundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(X(1:scarseInd:end,1:scarseInd:end)./1000,Y(1:scarseInd:end,1:scarseInd:end)./1000,tmpIZ(1:scarseInd:end,1:scarseInd:end),'.r')

        xlabel('X (km)')
        ylabel('Y (km)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        if ~legendHit
            legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
            legendHit = 1;
        end
        title([figPre,'IDW Interpolation'])
        view(20,20);
    end

    if GI.inputs.interpSurfaces.Kriging
        [tmpIZ,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Kriging);
        subplot(2,2,3)
        hold on
        plot3(edXYZ(:,1)./1000,edXYZ(:,2)./1000,edXYZ(:,3),'.k')
        plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,GP.boundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(X(1:scarseInd:end,1:scarseInd:end)./1000,Y(1:scarseInd:end,1:scarseInd:end)./1000,tmpIZ(1:scarseInd:end,1:scarseInd:end),'.r')

        xlabel('X (km)')
        ylabel('Y (km)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        if ~legendHit
            legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
        end
        title([figPre,'Kriging Interpolation'])
        view(20,20);
    end
    savePlot(GI,gcf,'Basal_Surface')

catch
    warning('Could not plot Basal Surfaces')
end

%% Edifice Profiles - Non-Exaggerated
try
    [ii_IDL,jj_IDL] = find(Z == max(Z(:)));
    xZprof = Z(ii_IDL,:);
    yZprof = Z(:,jj_IDL);

    figure('Name','Scaled Profiles','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    subplot(2,1,1)
    hold on
    p1 = plot(x./1000,xZprof./1000,'-k','linewidth',2);
    p2 = plot([x(1),x(end)]./1000,[1,1]*GP.Summit_Contour./1000,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GI.inputs.interpSurfaces.Natural
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Natural);
        pp = plot(x./1000,tmpZI(ii_IDL,:)./1000,'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GI.inputs.interpSurfaces.IDW
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.IDW);
        pp = plot(x./1000,tmpZI(ii_IDL,:)./1000,'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GI.inputs.interpSurfaces.Kriging
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Kriging);
        pp = plot(x./1000,tmpZI(ii_IDL,:)./1000,'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kriging Interp.'}];
    end
    axis tight
    axis image
    i1 = find(~isnan(xZprof),1);
    i2 = find(~isnan(xZprof),1,'last');
    xlim([x(i1),x(i2)]./1000);

    box on
    xlabel('X (km)')
    ylabel('Z (km)')
    setAxes(gca,fs)
    title([figPre,'E-W Profile Through Highest Peak'])
    legend(allP,allPt,'location','best')

    subplot(2,1,2)
    hold on
    p1 = plot(y./1000,yZprof./1000,'-k','linewidth',2);
    p2 = plot([y(1),y(end)]./1000,[1,1]*GP.Summit_Contour./1000,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GI.inputs.interpSurfaces.Natural
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Natural);
        pp = plot(y./1000,tmpZI(:,jj_IDL)./1000,'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GI.inputs.interpSurfaces.IDW
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.IDW);
        pp = plot(y./1000,tmpZI(:,jj_IDL)./1000,'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GI.inputs.interpSurfaces.Kriging
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Kriging);
        pp = plot(y./1000,tmpZI(:,jj_IDL)./1000,'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kriging Interp.'}];
    end
    axis tight
    axis image
    i1 = find(~isnan(yZprof),1);
    i2 = find(~isnan(yZprof),1,'last');
    if y(i1) < y(i2)
        xlim([y(i1),y(i2)]./1000);
    else
        xlim([y(i2),y(i1)]./1000);
    end

    box on
    xlabel('Y (km)')
    ylabel('Z (km)')
    setAxes(gca,fs)
    title([figPre,'N-S Profile Through Highest Peak'])

    savePlot(GI,gcf,'Scaled_Profiles')
catch
    warning('Could not plot Non-Exaggerated Profiles')
end

%% Edifice Profiles - Exaggerated
try
    figure('Name','E-W Non-Scaled Profile','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    p1 = plot(x./1000,xZprof,'-k','linewidth',2);
    p2 = plot([x(1),x(end)]./1000,[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GI.inputs.interpSurfaces.Natural
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Natural);
        pp = plot(x./1000,tmpZI(ii_IDL,:),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GI.inputs.interpSurfaces.IDW
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.IDW);
        pp = plot(x./1000,tmpZI(ii_IDL,:),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GI.inputs.interpSurfaces.Kriging
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Kriging);
        pp = plot(x./1000,tmpZI(ii_IDL,:),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kriging Interp.'}];
    end

    axis tight

    i1 = find(~isnan(xZprof),1);
    i2 = find(~isnan(xZprof),1,'last');
    xlim([x(i1),x(i2)]./1000);

    box on
    xlabel('X (km)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title([figPre,'E-W Profile Through Highest Peak'])
    legend(allP,allPt)

    savePlot(GI,gcf,'NtS_EW_Profile')

    figure('Name','N-S Non-Scaled Profile','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    p1 = plot(y./1000,yZprof,'-k','linewidth',2);
    p2 = plot([y(1),y(end)]./1000,[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GI.inputs.interpSurfaces.Natural
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Natural);
        pp = plot(y./1000,tmpZI(:,jj_IDL),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GI.inputs.interpSurfaces.IDW
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.IDW);
        pp = plot(y./1000,tmpZI(:,jj_IDL),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GI.inputs.interpSurfaces.Kriging
        [tmpZI,~,~] = GRIDobj2mat(GP.Basal_Surfaces.Kriging);
        pp = plot(y./1000,tmpZI(:,jj_IDL),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kriging Interp.'}];
    end
    axis tight

    i1 = find(~isnan(yZprof),1);
    i2 = find(~isnan(yZprof),1,'last');
    if y(i1) < y(i2)
        xlim([y(i1),y(i2)]./1000);
    else
        xlim([y(i2),y(i1)]./1000);
    end

    box on
    xlabel('Y (km)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title([figPre,'N-S Profile Through Highest Peak'])
    legend(allP,allPt)

    savePlot(GI,gcf,'NtS_NS_Profile')

catch
    warning('Could not plot Exaggerated Profiles')
end

%% 3D Edifice Eroded Topography
try
    az = -37;
    el = 30;
    ErodeGrid = SiP.ReconstructedTopo.Convex_Hull_Interpolated_Surface - Z;
    ReconInter = scatteredInterpolant(X(:),Y(:),SiP.ReconstructedTopo.Convex_Hull_Interpolated_Surface(:));
    R_Bound_Z = ReconInter(GP.boundaryXYZ(:,1),GP.boundaryXYZ(:,2));
    R_Summit_Z = ReconInter(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2));
    R_LFlank_Z = ReconInter(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2));

    figure('Name','3D Eroded Topography','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    surf(X./1000,Y./1000,SiP.ReconstructedTopo.Convex_Hull_Interpolated_Surface,ErodeGrid); shading flat; view(az,el); colormap(bluewhitered(255)); camlight

    p1 = plot3(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,R_Bound_Z,'-k','linewidth',2);
    p2 = plot3(GP.SummitXYZ(:,1)./1000,GP.SummitXYZ(:,2)./1000,R_Summit_Z,'-r','linewidth',2);
    p4 = plot3(GP.Lower_Flank_XY(:,1)./1000,GP.Lower_Flank_XY(:,2)./1000,R_LFlank_Z,'--k','linewidth',2);
    title([figPre,'3D Reconstructed Surface'])

    legend([p1,p2,p4],{'Edifice Boundary','Summit Boundary','Lower Flank Boundary'});
    box on

    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Eroded Topography (m)')
    setAxes(cb,0)
    grid on
    axis tight

    savePlot(GI,gcf,'3D_ErodedTopography')
catch
    warning('Could not plot 3D Eroded Topography')
end

%% Edifice Eroded Topography Maps
try
    figure('Name','Eroded Topography Maps','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    subplot(2,2,1)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1);
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k','linewidth',1);
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title([figPre,sprintf('Current\nTopography')])
    box on

    subplot(2,2,2)
    hold on
    DEMR = GRIDobj(X,Y,SiP.ReconstructedTopo.Convex_Hull_Interpolated_Surface);
    imageschs(DEMR,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1);
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k','linewidth',1);
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title([figPre,sprintf('Reconstructed\nTopography')])
    box on

    subplot(2,2,[3,4])
    hold on
    DEME = GRIDobj(X,Y,SiP.ReconstructedTopo.Convex_Hull_Interpolated_Surface-Z);
    figure
    pcolor(DEME.Z);shading flat;
    cc = colormap(bluewhitered(255));
    close

    imageschs(GP.DEM,DEME,'colormap',cc,'tickstokm',1);
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k','linewidth',1);
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)

    cb = colorbar;
    ylabel(cb,'Eroded Topography (m)')
    setAxes(cb,0)
    title([figPre,sprintf('Eroded\nTopography')])
    box on

    savePlot(GI,gcf,'ErodedTopography')
catch
    warning('Could not plot Eroded Topography')
end

%% Edifice Contours
try
    figure('Name','Contours','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on

    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1)
    xlabel('X (km)')
    ylabel('Y (km)')

    [~,p1] = contour(X./1000,Y./1000,Z,ShP.Contour.Elevations,'-k','linewidth',2);
    r = 0:.01:2*pi;
    p3 = [];
    for i = 1:length(ShP.Contour.Ellipses)
        if ~isempty(ShP.Contour.Ellipses{i})
            CE = ShP.Contour.Ellipses{i};
            xFit = CE.longAxis*cos(r);
            yFit = CE.shortAxis*sin(r);
            rotMat = [cosd(-CE.mathPhi),-sind(-CE.mathPhi);sind(-CE.mathPhi),cosd(-CE.mathPhi)];
            xyRot = [xFit(:),yFit(:)]*rotMat;
            xFit = xyRot(:,1) + CE.x0;
            yFit = xyRot(:,2) + CE.y0;

            if ShP.Contour.Elevations(i) >= GP.Lower_Flank_Contour
                p2 = plot(xFit./1000,yFit./1000,'-r','linewidth',2);
            else
                p3 = plot(xFit./1000,yFit./1000,'--r','linewidth',2);
            end
        end
    end

    if isempty(p3)
        legend([p1,p2],'Elevation Contour','Best-Fitting Ellipse')
    else
        legend([p1,p2,p3],'Elevation Contour','Best-Fitting Ellipse','Open Contour Best-Fitting Ellipse')
    end
    title([figPre,'Analyzed Contours'])

    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)

    savePlot(GI,gcf,'Contours')

catch
    warning('Could not plot Contours')
end

%% Edifice Contour Statistics
try
    lw1 = 2;
    ms = 8;

    figure('Name','Contour Statistics','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    subplot(1,5,1)
    hold on
    for i = 1:size(SP.Contour.Elevations,1)
        p1 = plot(([SP.Contour.Min(i),SP.Contour.Max(i)]),...
            [1,1]*SP.Contour.Elevations(i),'-k','linewidth',1);
    end
    p2 = plot((SP.Contour.Mean),SP.Contour.Elevations,...
        '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
    p3 = plot((SP.Contour.Median),SP.Contour.Elevations(:,1),...
        '-bo','linewidth',1,'markerfacecolor','b','markeredgecolor','k','markersize',ms);
    xlabel('Slope (^o)')
    ylabel('Elevation (m)')
    setAxes(gca,fs)
    box on
    axis tight
    xx = xlim;
    p4 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p5 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour.Elevations),max(ShP.Contour.Elevations)])
    legend([p1,p2,p3,p4,p5],{'Slope Range','Mean Slope','Median Slope',...
        'Summit Contour','Minimum Closed Contour'},'location','southeast')
    title([figPre,'Contour Slopes'])
    
    subplot(1,5,2)
    hold on
    p3 = plot(ShP.Contour.Ellipticity_Indices.MaxDiameter,ShP.Contour.Elevations,...
        '-ko','markerfacecolor','k','markersize',ms);
    p4 = plot(ShP.Contour.Ellipticity_Indices.BFEllipse,ShP.Contour.Elevations,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k');
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour.Elevations),max(ShP.Contour.Elevations)])
    legend([p1,p2,p3,p4],{'Summit Contour','Lowest Closed Contour','Maximum Diameter','Best-Fitting Ellipse'},'location','southeast')
    xlabel('Ellipticity Index')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title([figPre,sprintf('Ellipticity Indexes')])
    
    subplot(1,5,3)
    hold on
    p3 = plot(ShP.Contour.Irregularity_Indices.MaxDiameter,ShP.Contour.Elevations,...
        '-ko','markerfacecolor','k','markersize',ms);
    p4 = plot(ShP.Contour.Irregularity_Indices.BFEllipse,ShP.Contour.Elevations,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k');
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour.Elevations),max(ShP.Contour.Elevations)])
%     legend([p1,p2,p3,p4],{'Summit Contour','Lowest Closed Contour','Maximum Diameter','Best-Fitting Ellipse'},'location','southeast')
    xlabel('Irregularity Index')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title([figPre,sprintf('Irregularity Indexes')])
    
    subplot(1,5,4)
    hold on
    plot(ShP.Contour.Axis_Ellipticity,ShP.Contour.Elevations,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k')
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour.Elevations),max(ShP.Contour.Elevations)])
%     legend([p1,p2],{'Summit Contour','Lowest Closed Contour'},'location','southeast')
    xlabel('Ellipse Ellipticity')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title([figPre,'Ellipse Ellipticity'])

    subplot(1,5,5)
    hold on
    plot(TP.SlopeVariance.Elevation.Values(:,5),TP.SlopeVariance.Elevation.Values(:,1),'-ok','linewidth',2,'markerfacecolor','k')
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    setAxes(gca,fs)
    xlabel('Slope Variance')
    box on
    axis tight
    set(gca,'ytick',[])
    ylim([min(TP.SlopeVariance.Elevation.Values(:,1)),max(TP.SlopeVariance.Elevation.Values(:,1))])
    title({[figPre,'Slope Variance'];sprintf('Total Slope Variance = %.2f',TP.SlopeVariance.Total)})

    savePlot(GI,gcf,'Edifice_Contour_Stat')
catch
    warning('Could not plot Contour Statistics')
end

%% Contour Roughness
try
    lw1 = 2;
    ms = 8;
    pCutoff = 5;

    figure('Name','Contour Roughness Statistics','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    if length(TP.Roughness) > pCutoff
        warning('Number of roughness windows exceeds figure panels. Only plotting %d.',pCutoff)
    end

    numP = min([pCutoff,length(TP.Roughness)]);
    for i = 1:numP
        subplot(1,pCutoff,i)
        hold on
        for j = 1:length(RP.Contour.Elevations)
            p1 = plot(([RP.Contour.Min(j,i),RP.Contour.Max(j,i)]),...
            [1,1]*RP.Contour.Elevations(j),'-k','linewidth',1);
        end

         p2 = plot((RP.Contour.Mean(:,i)),RP.Contour.Elevations,...
            '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
        p3 = plot((RP.Contour.Median(:,i)),RP.Contour.Elevations(:,1),...
            '-bo','linewidth',1,'markerfacecolor','b','markeredgecolor','k','markersize',ms);
        xlabel('Roughness')

        
        setAxes(gca,fs)
        box on
        axis tight
        xx = xlim;
        p4 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
        p5 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
        xlim(xx)
        ylim([min(RP.Contour.Elevations),max(RP.Contour.Elevations)])
        
        if i == 1
            ylabel('Elevation (m)')
            legend([p1,p2,p3,p4,p5],{'Roughness Range','Mean Roughness','Median Roughness',...
                'Summit Contour','Minimum Closed Contour'},'location','southeast')
        else
            set(gca,'ytick',[])
        end
        
        title([figPre,sprintf('Contour Roughness\n%.2f m Window Size',TP.Roughness(i).TrueWindowRes)])
    end

    savePlot(GI,gcf,'Contour_Roughness')
catch
    warning('Could not plot Contour Roughness')
end

%% Contour Windowed Slope Variance
try
    lw1 = 2;
    ms = 8;
    pCutoff = 5;

    figure('Name','Contour Windowed Slope Variance Statistics','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    if length(TP.SlopeVariance.Windows) > pCutoff
        warning('Number of slope variance windows exceeds figure panels. Only plotting %d.',pCutoff)
    end

    numP = min([pCutoff,length(TP.SlopeVariance.Windows)]);
    for i = 1:numP
        subplot(1,pCutoff,i)
        hold on
        for j = 1:length(WsvP.Contour.Elevations)
            p1 = plot(([WsvP.Contour.Min(j,i),WsvP.Contour.Max(j,i)]),...
            [1,1]*WsvP.Contour.Elevations(j),'-k','linewidth',1);
        end

         p2 = plot((WsvP.Contour.Mean(:,i)),WsvP.Contour.Elevations,...
            '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
        p3 = plot((WsvP.Contour.Median(:,i)),WsvP.Contour.Elevations(:,1),...
            '-bo','linewidth',1,'markerfacecolor','b','markeredgecolor','k','markersize',ms);
        xlabel('Windowed Slope Variance')

        
        setAxes(gca,fs)
        box on
        axis tight
        xx = xlim;
        p4 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
        p5 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
        xlim(xx)
        ylim([min(WsvP.Contour.Elevations),max(WsvP.Contour.Elevations)])
        
        if i == 1
            ylabel('Elevation (m)')
            legend([p1,p2,p3,p4,p5],{'Win. Slp. Var. Range','Mean Win. Slp. Var.','Median Win. Slp. Var.',...
                'Summit Contour','Minimum Closed Contour'},'location','southeast')
        else
            set(gca,'ytick',[])
        end
        
        title([figPre,sprintf('Contour Win. Slp. Var.\n%.2f m Window Size',TP.Roughness(i).TrueWindowRes)])
    end

    savePlot(GI,gcf,'Contour_Windowed_Slope_Variance')
catch
    warning('Could not plot Contour Windowed Slope Variance')
end

%% Peak Parameters
try
    ms1 = 14;
    ms2 = 12;

    figure('Name','Peak Parameters','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
    hold on
    imageschs(GP.DEM,[],'colormap',demcmap(GP.DEM.Z),'tickstokm',1)
    plot(GP.boundaryXYZ(:,1)./1000,GP.boundaryXYZ(:,2)./1000,'-k','linewidth',2);
    plot(GP.SummitXYZ(:,1)./1000,GP.SummitXYZ(:,2)./1000,'-r','linewidth',2);
    plot(GP.Lower_Flank_XY(:,1)./1000,GP.Lower_Flank_XY(:,2)./1000,'-b','linewidth',2);
    if ~isempty(GP.maskXY)
        for i = 1:length(GP.maskXY)
            fill(GP.maskXY{i}(:,1)./1000,GP.maskXY{i}(:,2)./1000,'w');
        end
    end
    if ~isempty(GP.craterXYZ)
        for i = 1:length(GP.craterXYZ)
            fill(GP.craterXYZ{i}(:,1)./1000,GP.craterXYZ{i}(:,2)./1000,'r','facealpha',.5);
        end
    end

    allP = [];
    allPT = {};
    % Summit Peaks
    if ~isempty(PP.Summit.Local.XYZ)
        pp = plot(PP.Summit.Local.XYZ(:,1)./1000,PP.Summit.Local.XYZ(:,2)./1000,'.r','markersize',ms1);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Summit Local Peak (%d)',PP.Summit.Local.Count)}];
    end
    if ~isempty(PP.Summit.Contour.XYZ)
        pp = plot(PP.Summit.Contour.XYZ(:,1)./1000,PP.Summit.Contour.XYZ(:,2)./1000,'ok','markerfacecolor','r','markersize',ms2);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Summit Contour Peak (%d)',PP.Summit.Contour.Count)}];
    end

    % Main Flank Peaks
    if ~isempty(PP.Main_Flank.Local.XYZ)
        pp = plot(PP.Main_Flank.Local.XYZ(:,1)./1000,PP.Main_Flank.Local.XYZ(:,2)./1000,'.b','markersize',ms1);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Main Flank Local Peak (%d)',PP.Main_Flank.Local.Count)}];
    end
    if ~isempty(PP.Main_Flank.Contour.XYZ)
        pp = plot(PP.Main_Flank.Contour.XYZ(:,1)./1000,PP.Main_Flank.Contour.XYZ(:,2)./1000,'ok','markerfacecolor','b','markersize',ms2);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Main Flank Contour Peak (%d)',PP.Main_Flank.Contour.Count)}];
    end

    % Lower Flank Peaks
    if ~isempty(PP.Lower_Flank.Local.XYZ)
        pp = plot(PP.Lower_Flank.Local.XYZ(:,1)./1000,PP.Lower_Flank.Local.XYZ(:,2)./1000,'.k','markersize',ms1);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Lower Flank Local Peak (%d)',PP.Lower_Flank.Local.Count)}];
    end
    if ~isempty(PP.Lower_Flank.Contour.XYZ)
        pp = plot(PP.Lower_Flank.Contour.XYZ(:,1)./1000,PP.Lower_Flank.Contour.XYZ(:,2)./1000,'ok','markerfacecolor','k','markersize',ms2);
        allP = [allP;pp];
        allPT = [allPT;{sprintf('Lower Flank Contour Peak (%d)',PP.Lower_Flank.Contour.Count)}];
    end
    
    xlabel('X (km)')
    ylabel('Y (km)')
    setAxes(gca,fs)
    cb = colorbar;
    setAxes(cb,[]);
    title({[figPre,'Peak Count'];...
        sprintf('(Total Local Peaks: %d)',PP.Full_Edifice.Local.Count);...
        sprintf('(Total Contour Peaks: %d)',PP.Full_Edifice.Contour.Count)})
    legend(allP,allPT,'Location','southwestoutside');
    savePlot(GI,gcf,'Peak_Parameters')
catch 
    warning('Could not plot Peak Parameters')
end

%% Crater Contour Stats
try
    if ~isempty(res.CraterParams)
        for i = 1:length(res.CraterParams.Slope.Contour)
            cc = res.CraterParams.Slope.Contour{i};
            
            figure('Name',sprintf('Crater %d Contour Statistics',i),'units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
            subplot(1,3,2)
            hold on
            for j = 1:size(cc.Elevations,1)
                p1 = plot([cc.Min(j),cc.Max(j)],...
                    [1,1]*cc.Elevations(j),...
                    '-k','linewidth',1);
            end
            p2 = plot(cc.Mean,cc.Elevations,...
                '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
            p3 = plot(cc.Median,cc.Elevations,...
                '-bo','markerfacecolor','b','markeredgecolor','k','markersize',ms);
            
            xlabel('Slope (^o)')
            ylabel('Elevation (m)')
            setAxes(gca,fs)
            legend([p1,p2,p3],{'Slope Range','Mean Slope','Median Slope'},'location','northeastoutside')
            title([figPre,sprintf('Crater %d Slopes',i)])
            axis tight
            box on
            
            savePlot(GI,gcf,sprintf('Crater_%d_Contour_Stats',i))
        end
    end
catch
    warning('Could not plot Crater Statistics')
end
    
%% Crater Surface Points
try
    scarseInd = 1;
    tX = X(1:scarseInd:end,1:scarseInd:end);
    tY = Y(1:scarseInd:end,1:scarseInd:end);
    tZ = Z(1:scarseInd:end,1:scarseInd:end);
    pointsXYZ = [tX(:),tY(:),tZ(:)];
    pointsXYZ(isnan(pointsXYZ(:,3)),:) = [];
    
    for i = 1:length(GP.craterXYZ)
        cxyz = GP.craterXYZ{i};
        tmpPoints = pointsXYZ;
        pp = inpolygon(tmpPoints(:,1),tmpPoints(:,2),cxyz(:,1),cxyz(:,2));
        tmpPoints(pp==0,:) = [];
        
        figure('Name','Crater Surface','units','normalized','outerposition',[0 0 1 1],'visible',GI.inputs.visPlots)
        
        legendHit = 0;
        if GI.inputs.interpSurfaces.Natural
            [tmpZI,~,~] = GRIDobj2mat(res.CraterParams.Crater_Surfaces{i}.Natural);
            subplot(2,2,1)
            hold on
            plot3(tmpPoints(:,1)./1000,tmpPoints(:,2)./1000,tmpPoints(:,3),'.k')
            plot3(cxyz(:,1)./1000,cxyz(:,2)./1000,cxyz(:,3),'ko','markerfacecolor','b')
            plot3(X./1000,Y./1000,tmpZI,'.r')

            xlabel('X (km)')
            ylabel('Y (km)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            legend('Crater Points','Boundary Points','Crater Surface Points')
            legendHit = 1;
            title([figPre,'Natural Interpolation'])
            view(20,20);
        end

        if GI.inputs.interpSurfaces.IDW
            [tmpZI,~,~] = GRIDobj2mat(res.CraterParams.Crater_Surfaces{i}.IDW);
            subplot(2,2,2)
            hold on
            plot3(tmpPoints(:,1)./1000,tmpPoints(:,2)./1000,tmpPoints(:,3),'.k')
            plot3(cxyz(:,1)./1000,cxyz(:,2)./1000,cxyz(:,3),'ko','markerfacecolor','b')
            plot3(X./1000,Y./1000,tmpZI,'.r')
            
            xlabel('X (km)')
            ylabel('Y (km)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            if ~legendHit
                legend('Crater Points','Boundary Points','Crater Surface Points')
                legendHit = 1;
            end
            title([figPre,'IDW Interpolation'])
            view(20,20);
        end

        if GI.inputs.interpSurfaces.Kriging
            [tmpZI,~,~] = GRIDobj2mat(res.CraterParams.Crater_Surfaces{i}.Kriging);
            subplot(2,2,3)
            hold on
            plot3(tmpPoints(:,1)./1000,tmpPoints(:,2)./1000,tmpPoints(:,3),'.k')
            plot3(cxyz(:,1)./1000,cxyz(:,2)./1000,cxyz(:,3),'ko','markerfacecolor','b')
            plot3(X./1000,Y./1000,tmpZI,'.r')
            
            xlabel('X (km)')
            ylabel('Y (km)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            if ~legendHit
                legend('Crater Points','Boundary Points','Crater Surface Points')
            end
            title([figPre,'Kriging Interpolation'])
            view(20,20);
        end

        savePlot(GI,gcf,sprintf('Crater %d_Surface',i))
    end 
catch
    warning('Could not plot Crater Surface')
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
        % saveas(h,[GP.inputs.saveFigFolder,GP.inputs.figPrefix,name,'.png']);
    end
end