function [Total_Chi_Class,Individual_Chi_Class,zCutoff] = DrainageVolc_Calculate_Chi(DEM0,DEM,DB,boundaryXY,craterXY,channelThreshold,topN_basinIDs,allMN0,plotResults,saveFigFolder,figPrefix,runParallel,zCutoff,removeUpperBasins,verbose)
% Name: DrainageVolc_Calculate_Chi
% Author: Daniel O'Hara
% Date: 06/09/2024 (mm/dd/yyyy)
% Description: Script to calculate Chi values. Unless already supplied, 
%   this algorihtm will first calculate an M/N value using optimization of
%   all basins, then use that value to determine Chi for all basins
%   together. Afterwards, individual M/N values will be deteremined and chi
%   calculated for each basin.
%
% Input:
%   DEM0: Initial GRIDobj of DEM.
%   DEM: Cut GRIDobj of DEM.
%   DB: GRIDobj of drainage basins.
%   boundaryXY: Array of boundary x-y coordinates.
%   craterXY: Cell array of crater x-y coordinates.
%   channelThreshold: Drainage area threshold for channelization.
%   topN_basinIDs: Basin IDs of the most important basins.
%   allMN0: Supplied M/N value. If NaN, M/N will be determined using top
%       basins.
%   plotResults: Flag to plot the results.
%   saveFigFolder: Folder path to save the figure.
%   figPrefix: Prefix to label the figure.
%   runParallel: Flag for whether to use parallelization.
%   chi_Zcutoff = Cutoff elevation for \Chi analysis. If it is set to 
%       NaN, the elevation is set to the highest boundary point. If set
%       to an imaginary number > 0 the elevation is set to that amount
%       of relief above the lowest edifice elevation (e.g., 50i will
%       set the elevation to 50 m above the lowest point). If set to an
%       imaginary number between -1 and 0, the cutoff is determined as
%       a percent of the total relief above the lowest point.
%   chi_removeUpperBasins = Flag for whether basins that have outlets 
%       above the elevation cutoff should be removed from the \Chi 
%       analysis.
%   verbose: Flag for whether updates should be outputted.
%   
% Output:
%   Total_Chi_Class: Structure that contains the cumulative basin chi
%       analysis. Contains:
%         Chi_S: STREAMobj of all edifice channel networks.
%         Chi_S_TopN: STREAMobj of only channel networks associated with
%           the top basins.
%         MN: Best-fitting M/N value (if not supplied).
%         Chi: List of Chi values for the Chi_S network.
%         Maximum_Chi: GRIDobj of maximum Chi values for each basin.
%         Upstream_Chi: GRIDobj of Chi values projected upstream to their 
%           associated divides. 
%   Inividual_Chi_Class: Structure that contains individual basin M/N and
%       Chi values. Contains:
%         Chi_S: Class of individual basin IDs and stream networks.
%         Chi: Class of individual basin IDs, M/N values, and Chi values.
%         Maximum_Chi: GRIDobj of maximum Chi values for each basin.
%         Upstream_Chi: GRIDobj of Chi values projected upstream to their 
%           associated divides. 
%   zCutoff: Chosen cutoff elevation.

%% Get minimum elevation for chi
if verbose > 0
    disp('      Finding minimum elevation for Chi...')
end

if isnan(zCutoff)
    [tmpZ,tmpX,tmpY] = GRIDobj2mat(DEM0);
    [tmpX,tmpY] = meshgrid(tmpX,tmpY);
    boundaryZ = interp2(tmpX,tmpY,tmpZ,boundaryXY(:,1),boundaryXY(:,2));
    zCutoff = max(boundaryZ);
elseif ~isreal(zCutoff)
    zCutoff = imag(zCutoff);
    if zCutoff > 0
        zCutoff = nanmin(DEM.Z(:)) + zCutoff;
    elseif zCutoff >= -1 && zCutoff <= 0
        zCutoff = nanmin(DEM.Z(:)) + range(DEM.Z(:))*zCutoff;
    end
end

%% Create template DEM
if verbose > 0
    disp('      Creating template GRIDobj and STREAMobj...')
end
templateChiDEM = DEM;
[ZG,XG,YG] = GRIDobj2mat(DEM);
[XG,YG] = meshgrid(XG,YG);

% Remove elevations less than the cutoff.
templateChiDEM.Z(templateChiDEM.Z<zCutoff) = NaN;

% Remove the crater
for i = 1:length(craterXY)
    p = inpolygon(XG,YG,craterXY{i}(:,1),craterXY{i}(:,2));
    templateChiDEM.Z(p) = NaN;
end

% Remove basins with outlets higher than the cutoff
[DBg,~,~] = GRIDobj2mat(DB);
DBg(isnan(ZG)) = NaN;
uniDB = unique(DBg(:));
uniDB(uniDB==0) = NaN;
uniDB(isnan(uniDB)) = [];

if removeUpperBasins
    for i = 1:length(uniDB)
        minZ = templateChiDEM.Z(DB.Z==uniDB(i));
        if min(minZ(:)) > zCutoff + range(DEM.Z(:))*.01
            templateChiDEM.Z(DB.Z==uniDB(i)) = NaN;
        end
    end
end

templateChiDEM_allBasins = templateChiDEM;
FD_allBasins = FLOWobj(templateChiDEM_allBasins);
A_allBasins = flowacc(FD_allBasins);
totalChiS = STREAMobj(FD_allBasins,'minarea',channelThreshold,'unit','mapunits');

% Remove non top basins
for i = 1:length(uniDB)
    if sum(topN_basinIDs==uniDB(i))
        continue;
    end

    templateChiDEM.Z(DB.Z==uniDB(i)) = NaN;
end

%% Calculate MN and chi for the entire edifice
if verbose > 0
    disp('      Calculating Chi for the entire edifice...')
end
FD = FLOWobj(templateChiDEM,'preprocess','none');
A = flowacc(FD);
totalChiS_topN = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');
allMN = allMN0;

if ~isempty(totalChiS_topN.ix)
    try
        if isnan(allMN)
            if ~plotResults
                try
                    [allMN,~,~,~] = mnoptimvar(totalChiS_topN,templateChiDEM,A,'varfun',@robustcov,'plot',false);
                catch
                    try
                        [allMN,~,~,~] = mnoptimvar(totalChiS_topN,templateChiDEM,A,'plot',false);
                    catch
                        set(0,'DefaultFigureVisible','off')
                        try
                            [allMN,~] = mnoptim(totalChiS_topN,templateChiDEM,A,'mnrange',[.01,4],'UseParallel',runParallel);close; close;
                            allMN = allMN.mn;
                        catch
                            [allMN,~] = mnoptim(totalChiS_topN,templateChiDEM,A,'mnrange',[.01,4],'UseParallel',runParallel,'crossval',false);close; close
                            allMN = allMN.mn;
                        end
                        set(0,'DefaultFigureVisible','on')
                    end
                end
            else
                try
                    [allMN,~,~,~] = mnoptimvar(totalChiS_topN,templateChiDEM,A,'varfun',@robustcov);
                catch
                    try
                        [allMN,~,~,~] = mnoptimvar(totalChiS_topN,templateChiDEM,A);
                    catch
                        try
                            [allMN,~] = mnoptim(totalChiS_topN,templateChiDEM,A,'mnrange',[.01,4],'UseParallel',runParallel);close;
                            allMN = allMN.mn;
                        catch
                            [allMN,~] = mnoptim(totalChiS_topN,templateChiDEM,A,'mnrange',[.01,4],'UseParallel',runParallel,'crossval',false);close
                            allMN = allMN.mn;
                        end
                    end
                end
                if ~isempty(saveFigFolder)
                    saveas(gcf,[saveFigFolder,figPrefix,'Total_MN_Optimization.fig']);
                    saveas(gcf,[saveFigFolder,figPrefix,'Total_MN_Optimization.png']);
                end
                close;
            end
        end

        Total_Chi = chitransform(totalChiS,A_allBasins,'a0',channelThreshold,'mn',allMN);
        ChiG = STREAMobj2GRIDobj(totalChiS,Total_Chi);
        Total_DB_MaxChi = DB;
        Total_DB_MaxChi.Z(:) = NaN;
        for i = 1:length(uniDB)
                tmpDB = DB.Z == uniDB(i);
                maxVal = max(max(ChiG.Z(tmpDB==1)));
                
                Total_DB_MaxChi.Z(tmpDB==1) = maxVal;
        end

        Total_Upstream_Chi = DrainageVolc_Project_Chi_Upstream(DEM,DB,totalChiS,Total_Chi);
    catch
        allMN = allMN0;
        Total_Chi = NaN;
        Total_DB_MaxChi = NaN;
        Total_Upstream_Chi = NaN;
    end
else
    allMN = allMN0;
    Total_Chi = NaN;
    Total_DB_MaxChi = NaN;
    Total_Upstream_Chi = NaN;
end

Total_Chi_Class.Chi_S = totalChiS;
Total_Chi_Class.Chi_S_TopN = totalChiS_topN;
Total_Chi_Class.MN = allMN;
Total_Chi_Class.Chi = Total_Chi;
Total_Chi_Class.Maximum_Chi = Total_DB_MaxChi;
Total_Chi_Class.Upstream_Chi = Total_Upstream_Chi;

%% Calcualte MN and Chi for individual basins
if verbose > 0
    disp('      Calculating Chi for individual basins...')
end
warning('off','all')
set(0,'DefaultFigureVisible','off')

Individual_ChiS = [];
Individual_Chi = [];
Individual_DB_MaxChi = DB;
Individual_DB_MaxChi.Z(:) = NaN;
Individual_Upstream_Chi = DB;
Individual_Upstream_Chi.Z(:) = NaN;

curPer = 0;
counter = 1;
for i = 1:length(uniDB)
    if verbose > 0 && i/length(uniDB)>=curPer
        disp(sprintf('         %d%% Complete (%d / %d)',round(curPer*100,0),i,length(uniDB)))
        curPer = curPer + .1;
    end

    tmpDEM = templateChiDEM_allBasins;
    tmpDEM.Z(DB.Z~=uniDB(i)) = NaN;
    if sum(~isnan(tmpDEM.Z)) == 0
        continue;
    end
    FD = FLOWobj(tmpDEM,'preprocess','none');
    A = flowacc(FD);
    indiChiS = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');

    if isempty(indiChiS.ix)
        continue;
    end

    try
        [indiMN,~,~,~] = mnoptimvar(indiChiS,tmpDEM,A,'varfun',@robustcov,'plot',false);
    catch
        try
            [indiMN,~,~,~] = mnoptimvar(indiChiS,tmpDEM,A,'plot',false);
        catch
            try
                [indiMN,~] = mnoptim(indiChiS,tmpDEM,A,'mnrange',[.01,4],'UseParallel',runParallel);close; close
                indiMN = indiMN.mn;
            catch
                try
                    [indiMN,~] = mnoptim(indiChiS,tmpDEM,A,'mnrange',[.01,4],'UseParallel',runParallel,'crossval',false); close; close
                    indiMN = indiMN.mn;
                catch
                    continue;
                end
            end
        end
    end

    indiChi = chitransform(indiChiS,A,'a0',channelThreshold,'mn',indiMN);

    Individual_ChiS(counter).BasinID = uniDB(i);
    Individual_ChiS(counter).S = indiChiS;
    Individual_Chi(counter).BasinID = uniDB(i);
    Individual_Chi(counter).MN = indiMN;
    Individual_Chi(counter).Chi = indiChi;

    Individual_DB_MaxChi.Z(DB.Z==uniDB(i)) = nanmax(indiChi);
    tmpUpstreamChi = DrainageVolc_Project_Chi_Upstream(DEM,DB,indiChiS,indiChi);
    Individual_Upstream_Chi.Z(DB.Z==uniDB(i)) = tmpUpstreamChi.Z(DB.Z==uniDB(i));

    counter = counter + 1;
end
set(0,'DefaultFigureVisible','on')

Individual_Chi_Class.Chi_S = Individual_ChiS;
Individual_Chi_Class.Chi = Individual_Chi;
Individual_Chi_Class.Maximum_Chi = Individual_DB_MaxChi;
Individual_Chi_Class.Upstream_Chi = Individual_Upstream_Chi;
warning('on','all')

end