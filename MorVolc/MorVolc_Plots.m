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
GP = res.GeneralParams;
OP = res.OrientParams;
SP = res.SlopeParams;
ShP = res.ShapeParams;
SiP = res.SizeParams;

useDEM = viridis(255);

fs = 15;

%% Hillshade DEM
    figure('Name','DEM Hillshade','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(GP.DEM,[],'colormap',[.8 .8 .8],'colorbar',false)
    p1 = plot(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),'-k','linewidth',2);
    p2 = plot(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2),'-r','linewidth',2);
    p4 = plot(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2),'-b','linewidth',2);
    if ~isempty(GP.MaskXY)
        for i = 1:length(GP.MaskXY)
            p3 = fill(GP.MaskXY{i}(:,1),GP.MaskXY{i}(:,2),'w');
        end
        allP = [p1,p2,p4,p3];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary','Mask Region'};
    else
        allP = [p1,p2,p4];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary'};
    end
    if ~isempty(GP.CraterXYZ)
        for i = 1:length(GP.CraterXYZ)
            pp = fill(GP.CraterXYZ{i}(:,1),GP.CraterXYZ{i}(:,2),'r','facealpha',.5);
        end
        allP = [allP,pp];
        allPT = [allPT,{'Crater Region'}];
    end
    
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    title('Edifice Hillshade')
    legend(allP,allPT);
    savePlot(GP,gcf,'Hillshade_DEM')
    
%% 3D DEM
    figure('Name','3D DEM','units','normalized','outerposition',[0 0 1 1])
    surf(GP.DEM,'exaggerate',2); colormap(demcmap(GP.DEM.Z)); camlight
    hold on
    p1 = plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),GP.BoundaryXYZ(:,3),'-k','linewidth',2);
    p2 = plot3(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2),GP.SummitXYZ(:,3),'-r','linewidth',2);
    p4 = plot3(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2),GP.Lower_Flank_XY(:,3),'-b','linewidth',2);
    title('3D Topography (2x Exaggeration)')
%     legend([p1,p2,p4],{'Edifice Boundary','Summit Boundary','Lower Flank Boundary'});
    
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    savePlot(GP,gcf,'3D_DEM')
    
%% 2D Slope
    SDEM = GRIDobj(GP.X,GP.Y,GP.S);
    cDEM = GRIDobj(GP.X,GP.Y,GP.Z);
    figure('Name','2D Slope','units','normalized','outerposition',[0 0 1 1])
    hold on
    imageschs(cDEM,SDEM,'colormap',useDEM);
    p1 = plot(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),'-k','linewidth',1);
    p2 = plot(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2),'-r','linewidth',2);
    p4 = plot(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2),'--k','linewidth',2);
    if ~isempty(GP.MaskXY)
        for i = 1:length(GP.MaskXY)
            p3 = fill(GP.MaskXY{i}(:,1),GP.MaskXY{i}(:,2),'w');
        end
        allP = [p1,p2,p4,p3];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary','Mask Region'};
    else
        allP = [p1,p2,p4];
        allPT = {'Edifice Boundary','Summit Boundary','Lower Flank Boundary'};
    end
    if ~isempty(GP.CraterXYZ)
        for i = 1:length(GP.CraterXYZ)
            pp = plot(GP.CraterXYZ{i}(:,1),GP.CraterXYZ{i}(:,2),'-w','linewidth',2);
        end
        allP = [allP,pp];
        allPT = [allPT,{'Crater Boundary'}];
    end
    
    box on
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    title('Edifice Slopes')
    legend(allP,allPT);
    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    savePlot(GP,gcf,'2D_Slope')

%% 3D Slope
    figure('Name','3D Slope','units','normalized','outerposition',[0 0 1 1])
    surf(GP.X,GP.Y,GP.Z,GP.S); shading flat; colormap(useDEM); camlight
    hold on
    p1 = plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),GP.BoundaryXYZ(:,3),'-k','linewidth',2);
    p2 = plot3(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2),GP.SummitXYZ(:,3),'-r','linewidth',2);
    p4 = plot3(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2),GP.Lower_Flank_XY(:,3),'--k','linewidth',2);
    title('3D Slope (2x Exaggeration)')
    axis(gca,'image')
    set(gca,'DataAspectRatio',[1,1,1/2]);
    
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Slope (^o)')
    setAxes(cb,0)
    savePlot(GP,gcf,'3D_Slope')
    
%% Basal Surface Points
    scarseInd = 5;
    tX = GP.X(1:scarseInd:end,1:scarseInd:end);
    tY = GP.Y(1:scarseInd:end,1:scarseInd:end);
    tZ = GP.Z(1:scarseInd:end,1:scarseInd:end);
    pointsXYZ = [tX(:),tY(:),tZ(:)];
    pointsXYZ(isnan(pointsXYZ(:,3)),:) = [];
    
    edRel = range(pointsXYZ(:,3));
    edXYZ = pointsXYZ;
    edXYZ(edXYZ(:,3)>min(edXYZ(:,3))+edRel*.4,:) = [];
    
    figure('Name','Basal Surface','units','normalized','outerposition',[0 0 1 1])
    scarseInd = scarseInd*4;
    legendHit = 0;
    if GP.inputs.interpSurfaces.Natural
        subplot(2,2,1)
        hold on
        plot3(edXYZ(:,1),edXYZ(:,2),edXYZ(:,3),'.k')
        plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),GP.BoundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(GP.Basal_Surface_XY(1:scarseInd:end,1),GP.Basal_Surface_XY(1:scarseInd:end,2),GP.Basal_Surface_Z.Natural(1:scarseInd:end),'.r')
    
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        title('Basal Surface')
        legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
        legendHit = 1;
        title('Natural Interpolation')
        view(20,20);
    end
    
    if GP.inputs.interpSurfaces.IDW
        subplot(2,2,2)
        hold on
        plot3(edXYZ(:,1),edXYZ(:,2),edXYZ(:,3),'.k')
        plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),GP.BoundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(GP.Basal_Surface_XY(1:scarseInd:end,1),GP.Basal_Surface_XY(1:scarseInd:end,2),GP.Basal_Surface_Z.IDW(1:scarseInd:end),'.r')
    
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        title('Basal Surface')
        if ~legendHit
            legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
            legendHit = 1;
        end
        title('IDW Interpolation')
        view(20,20);
    end
    
    if GP.inputs.interpSurfaces.Kringing
        subplot(2,2,3)
        hold on
        plot3(edXYZ(:,1),edXYZ(:,2),edXYZ(:,3),'.k')
        plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),GP.BoundaryXYZ(:,3),'ko','markerfacecolor','b')
        plot3(GP.Basal_Surface_XY(1:scarseInd:end,1),GP.Basal_Surface_XY(1:scarseInd:end,2),GP.Basal_Surface_Z.Kringing(1:scarseInd:end),'.r')
    
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        setAxes(gca,fs)
        box on
        title('Basal Surface')
        if ~legendHit
            legend('Lower Edifice Points','Boundary Points','Basal Surface Points')
        end
        title('Kringing Interpolation')
        view(20,20);
    end
    savePlot(GP,gcf,'Basal_Surface')
    
%% Edifice Profiles - Non-Exaggerated
    [ii_IDL,jj_IDL] = find(GP.Z == max(GP.Z(:)));
    xTmp = GP.X;
    yTmp = GP.Y;
    xTmp(isnan(GP.Z)) = NaN;
    yTmp(isnan(GP.Z)) = NaN;
    xProf = xTmp(ii_IDL,:);
    yProf = yTmp(:,jj_IDL);
    
    xZprof = GP.Z(ii_IDL,:);
    yZprof = GP.Z(:,jj_IDL);
    
    xProfBPointsLog = GP.Basal_Surface_XY(:,2)==GP.Y(ii_IDL,jj_IDL);
    yProfBPointsLog = GP.Basal_Surface_XY(:,1)==GP.X(ii_IDL,jj_IDL);
    xProfBPoints = GP.Basal_Surface_XY(xProfBPointsLog,1);
    yProfBPoints = GP.Basal_Surface_XY(yProfBPointsLog,2);
    
    figure('Name','Scaled Profiles','units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1)
    hold on
    p1 = plot(xProf,xZprof,'-k','linewidth',2);
    p2 = plot([xProf(find(~isnan(xProf),1)),xProf(find(~isnan(xProf),1,'last'))],[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GP.inputs.interpSurfaces.Natural
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.Natural(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GP.inputs.interpSurfaces.IDW
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.IDW(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GP.inputs.interpSurfaces.Kringing
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.Kringing(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kringing Interp.'}];
    end
    axis tight
    axis image
    box on
    xlabel('X (m)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title('E-W Profile Through Highest Peak')
    legend(allP,allPt,'location','best')
    
    subplot(2,1,2)
    hold on
    p1 = plot(yProf,yZprof,'-k','linewidth',2);
    p2 = plot([yProf(find(~isnan(yProf),1)),yProf(find(~isnan(yProf),1,'last'))],[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GP.inputs.interpSurfaces.Natural
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.Natural(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GP.inputs.interpSurfaces.IDW
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.IDW(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GP.inputs.interpSurfaces.Kringing
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.Kringing(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kringing Interp.'}];
    end
    axis tight
    axis image
    box on
    xlabel('Y (m)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title('N-S Profile Through Highest Peak')
%     legend(allP,allPt)

    savePlot(GP,gcf,'Scaled_Profiles')
    
%% Edifice Profiles - Exaggerated
    figure('Name','E-W Non-Scaled Profile','units','normalized','outerposition',[0 0 1 1])
    hold on
    p1 = plot(xProf,xZprof,'-k','linewidth',2);
    p2 = plot([xProf(find(~isnan(xProf),1)),xProf(find(~isnan(xProf),1,'last'))],[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GP.inputs.interpSurfaces.Natural
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.Natural(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GP.inputs.interpSurfaces.IDW
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.IDW(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GP.inputs.interpSurfaces.Kringing
        tx = sortrows([xProfBPoints,GP.Basal_Surface_Z.Kringing(xProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kringing Interp.'}];
    end
    axis tight
    box on
    xlabel('X (m)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title('E-W Profile Through Highest Peak')
    legend(allP,allPt)
    
    savePlot(GP,gcf,'NtS_EW_Profile')
    
    
    figure('Name','N-S Non-Scaled Profile','units','normalized','outerposition',[0 0 1 1])
    hold on
    p1 = plot(yProf,yZprof,'-k','linewidth',2);
    p2 = plot([yProf(find(~isnan(yProf),1)),yProf(find(~isnan(yProf),1,'last'))],[1,1]*GP.Summit_Contour,'--k','linewidth',1);
    allP = [p1,p2];
    allPt = {'Topography','Summit Height'};
    if GP.inputs.interpSurfaces.Natural
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.Natural(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Natural Interp.'}];
    end
    if GP.inputs.interpSurfaces.IDW
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.IDW(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'-b','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'IDW Interp.'}];
    end
    if GP.inputs.interpSurfaces.Kringing
        tx = sortrows([yProfBPoints,GP.Basal_Surface_Z.Kringing(yProfBPointsLog)],1);
        pp = plot(tx(:,1),tx(:,2),'--r','linewidth',2);
        allP = [allP,pp];
        allPt = [allPt,{'Kringing Interp.'}];
    end
    axis tight
    box on
    xlabel('Y (m)')
    ylabel('Z (m)')
    setAxes(gca,fs)
    title('N-S Profile Through Highest Peak')
    legend(allP,allPt)
    
    savePlot(GP,gcf,'NtS_NS_Profile')
    
%% 3D Edifice Eroded Topography
    az = -37;
    el = 30;
    ErodeGrid = SiP.Convex_Hull_Interpolated_Surface - GP.Z;
    ReconInter = scatteredInterpolant(GP.X(:),GP.Y(:),SiP.Convex_Hull_Interpolated_Surface(:));
    R_Bound_Z = ReconInter(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2));
    R_Summit_Z = ReconInter(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2));
    R_LFlank_Z = ReconInter(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2));

    figure('Name','3D Eroded Topography','units','normalized','outerposition',[0 0 1 1])
    
    hold on
    surf(GP.X,GP.Y,SiP.Convex_Hull_Interpolated_Surface,ErodeGrid); shading flat; view(az,el); colormap(bluewhitered(255)); camlight
    
    p1 = plot3(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),R_Bound_Z,'-k','linewidth',2);
    p2 = plot3(GP.SummitXYZ(:,1),GP.SummitXYZ(:,2),R_Summit_Z,'-r','linewidth',2);
    p4 = plot3(GP.Lower_Flank_XY(:,1),GP.Lower_Flank_XY(:,2),R_LFlank_Z,'--k','linewidth',2);
    title('3D Reconstructed Surface')
    axis(gca,'image')
    set(gca,'DataAspectRatio',[1,1,1/2]);

    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Elevation (m)')
    setAxes(gca,fs)
    cb = colorbar;
    ylabel(cb,'Eroded Topography (m)')
    setAxes(cb,0)
    
    savePlot(GP,gcf,'3D_ErodedTopography')

%% 3D Edifice Eroded Topography
    figure('Name','Eroded Topography Maps','units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    hold on
    DEM0 = GRIDobj(GP.X,GP.Y,GP.Z);
    imageschs(DEM0,[],'colormap',demcmap(GP.DEM.Z));
    plot(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),'-k','linewidth',1);
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title(sprintf('Current\nTopography'))
    box on
    
    subplot(2,2,2)
    hold on
    DEMR = GRIDobj(GP.X,GP.Y,SiP.Convex_Hull_Interpolated_Surface);
    imageschs(DEMR,[],'colormap',demcmap(GP.DEM.Z));
    plot(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),'-k','linewidth',1);
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Elevation (m)')
    setAxes(cb,0)
    title(sprintf('Reconstructed\nTopography'))
    box on
    
    subplot(2,2,[3,4])
    hold on
    DEME = GRIDobj(GP.X,GP.Y,SiP.Convex_Hull_Interpolated_Surface-GP.Z);
    figure
    pcolor(DEME.Z);shading flat;
    cc = colormap(bluewhitered(255));
    close
    
    imageschs(DEM0,DEME,'colormap',cc);
    plot(GP.BoundaryXYZ(:,1),GP.BoundaryXYZ(:,2),'-k','linewidth',1);
    xlabel('X (m)')
    ylabel('Y (m)')
    setAxes(gca,fs)
    
    cb = colorbar;
    ylabel(cb,'Eroded Topography (m)')
    setAxes(cb,0)
    title(sprintf('Eroded\nTopography'))
    box on

    savePlot(GP,gcf,'ErodedTopography')

%% Edifice Contour Statistics
    lw1 = 2;
    ms = 8;

    figure('Name','Contour Statistics','units','normalized','outerposition',[0 0 1 1])
    subplot(1,4,1)
    hold on
    for i = 1:size(OP.Contour_Values,1)
        p1 = plot(([SP.Contour_Min_Slopes(i),SP.Contour_Max_Slopes(i)]),...
            [1,1]*SP.Contour_Values(i),'-k','linewidth',1);
    end
    p2 = plot((SP.Contour_Mean_Slopes),SP.Contour_Values,...
        '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
    p3 = plot((SP.Contour_Median_Slopes),SP.Contour_Values(:,1),...
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
    ylim([min(ShP.Contour_Values),max(ShP.Contour_Values)])
    legend([p1,p2,p3,p4,p5],{'Slope Range','Mean Slope','Median Slope',...
        'Summit Contour','Minimum Closed Contour'},'location','southeast')
    title('Contour Slopes')
    
    subplot(1,4,2)
    hold on
    p3 = plot(ShP.Contour_MaxDiam_Ellipticity_Indexes,ShP.Contour_Values,...
        '-ko','markerfacecolor','k','markersize',ms);
    p4 = plot(ShP.Contour_BFEllipse_Ellipticity_Indexes,ShP.Contour_Values,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k');
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour_Values),max(ShP.Contour_Values)])
    legend([p1,p2,p3,p4],{'Summit Contour','Lowest Closed Contour','Maximum Diameter','Best-Fitting Ellipse'},'location','southeast')
    xlabel('Ellipticity Index')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title(sprintf('Ellipticity Indexes'))
    
    subplot(1,4,3)
    hold on
    p3 = plot(ShP.Contour_MaxDiam_Irregularity_Indexes,ShP.Contour_Values,...
        '-ko','markerfacecolor','k','markersize',ms);
    p4 = plot(ShP.Contour_BFEllipse_Irregularity_Indexes,ShP.Contour_Values,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k');
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour_Values),max(ShP.Contour_Values)])
%     legend([p1,p2,p3,p4],{'Summit Contour','Lowest Closed Contour','Maximum Diameter','Best-Fitting Ellipse'},'location','southeast')
    xlabel('Irregularity Index')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title(sprintf('Irregularity Indexes'))
    
    subplot(1,4,4)
    hold on
    plot(ShP.Contour_Axis_Ellipticity,ShP.Contour_Values,...
        '-ro','markerfacecolor','r','markersize',ms,'markeredgecolor','k')
    box on
    axis tight
    xx = xlim;
    p1 = plot(xx,[1,1]*GP.Summit_Contour,'--r','linewidth',lw1);
    p2 = plot(xx,[1,1]*GP.Lower_Flank_Contour,'--b','linewidth',lw1);
    xlim(xx)
    ylim([min(ShP.Contour_Values),max(ShP.Contour_Values)])
%     legend([p1,p2],{'Summit Contour','Lowest Closed Contour'},'location','southeast')
    xlabel('Ellipse Ellipticity')
%     ylabel('Elevation (m)')
    set(gca,'ytick',[])
    setAxes(gca,fs)
    title('Ellipse Ellipticity')
    savePlot(GP,gcf,'Edifice_Contour_Stats')
    
%% Crater Contour Stats
    if ~isempty(res.CraterParams)
        for i = 1:length(res.CraterParams.Crater_Contour_Stats)
            cc = res.CraterParams.Crater_Contour_Stats{i};
            
            figure('Name',sprintf('Crater %d Contour Statistics',i),'units','normalized','outerposition',[0 0 1 1])
            subplot(1,3,2)
            hold on
            for j = 1:size(cc.Crater_Contours,1)
                p1 = plot([cc.Crater_Min_Slope(j),cc.Crater_Max_Slope(j)],...
                    [1,1]*cc.Crater_Contours(j),...
                    '-k','linewidth',1);
            end
            p2 = plot(cc.Crater_Mean_Slope,cc.Crater_Contours,...
                '-sr','linewidth',1,'markerfacecolor','r','markeredgecolor','k','markersize',ms);
            p3 = plot(cc.Crater_Median_Slope,cc.Crater_Contours,...
                '-bo','markerfacecolor','b','markeredgecolor','k','markersize',ms);
            
            xlabel('Slope (^o)')
            ylabel('Elevation (m)')
            setAxes(gca,fs)
            legend([p1,p2,p3],{'Slope Range','Mean Slope','Median Slope'},'location','northeastoutside')
            title(sprintf('Crater %d Slopes',i))
            axis tight
            box on
            
            savePlot(GP,gcf,sprintf('Crater_%d_Contour_Stats',i))
        end
    end
    
%% Crater Surface Points
    scarseInd = 1;
    tX = GP.X(1:scarseInd:end,1:scarseInd:end);
    tY = GP.Y(1:scarseInd:end,1:scarseInd:end);
    tZ = GP.Z(1:scarseInd:end,1:scarseInd:end);
    pointsXYZ = [tX(:),tY(:),tZ(:)];
    pointsXYZ(isnan(pointsXYZ(:,3)),:) = [];
    
    for i = 1:length(GP.CraterXYZ)
        cxyz = GP.CraterXYZ{i};
        tmpPoints = pointsXYZ;
        pp = inpolygon(tmpPoints(:,1),tmpPoints(:,2),cxyz(:,1),cxyz(:,2));
        tmpPoints(pp==0,:) = [];
        
        figure('Name','Crater Surface','units','normalized','outerposition',[0 0 1 1])
        
        legendHit = 0;
        if GP.inputs.interpSurfaces.Natural
            subplot(2,2,1)
            hold on
            plot3(tmpPoints(:,1),tmpPoints(:,2),tmpPoints(:,3),'.k')
            plot3(cxyz(:,1),cxyz(:,2),cxyz(:,3),'ko','markerfacecolor','b')
            plot3(res.CraterParams.Crater_Surface_XY{i}(:,1),res.CraterParams.Crater_Surface_XY{i}(:,2),res.CraterParams.Crater_Surface_Z{i}.Natural,'.r')

            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            legend('Crater Points','Boundary Points','Crater Surface Points')
            legendHit = 1;
            title('Natural Interpolation')
            view(20,20);
        end

        if GP.inputs.interpSurfaces.IDW
            subplot(2,2,2)
            hold on
            plot3(tmpPoints(:,1),tmpPoints(:,2),tmpPoints(:,3),'.k')
            plot3(cxyz(:,1),cxyz(:,2),cxyz(:,3),'ko','markerfacecolor','b')
            plot3(res.CraterParams.Crater_Surface_XY{i}(:,1),res.CraterParams.Crater_Surface_XY{i}(:,2),res.CraterParams.Crater_Surface_Z{i}.IDW,'.r')
            
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            if ~legendHit
                legend('Crater Points','Boundary Points','Crater Surface Points')
                legendHit = 1;
            end
            title('IDW Interpolation')
            view(20,20);
        end

        if GP.inputs.interpSurfaces.Kringing
            subplot(2,2,3)
            hold on
            plot3(tmpPoints(:,1),tmpPoints(:,2),tmpPoints(:,3),'.k')
            plot3(cxyz(:,1),cxyz(:,2),cxyz(:,3),'ko','markerfacecolor','b')
            plot3(res.CraterParams.Crater_Surface_XY{i}(:,1),res.CraterParams.Crater_Surface_XY{i}(:,2),res.CraterParams.Crater_Surface_Z{i}.Kringing,'.r')
            
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            setAxes(gca,fs)
            box on
            title('Basal Surface')
            if ~legendHit
                legend('Crater Points','Boundary Points','Crater Surface Points')
            end
            title('Kringing Interpolation')
            view(20,20);
        end

        savePlot(GP,gcf,sprintf('Crater %d_Surface',i))
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