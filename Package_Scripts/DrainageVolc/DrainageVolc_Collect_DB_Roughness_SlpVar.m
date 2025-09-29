function [db_rVals,db_svVals] = DrainageVolc_Collect_DB_Roughness_SlpVar(DB,rVals,svVals,verbose)
% Name: DrainageVolc_Collect_DB_Roughness_SlpVar
% Author: Daniel O'Hara
% Date: 06/07/2024 (mm/dd/yyyy)
% Description: Script to calculate mean basin roughnesses and slope
%   variances for a given set of windows 
%
% Input:
%   DB: GRIDobj of drainage basins.
%   rVals: Roughness structure. Given from the TopoParams class.
%   svVals: Slope Variance structure. Given from the TopoParams class.
%   verbose: Flag for outputting.
%
% Output:
%   db_rVals: Structure of roughness values, contains basin IDs and 
%       roughness windows, as well as mean, median, std, min, and max 
%       roughnesses. Statistics matrices have rows the same number as 
%       basins, and columns the same number as roughness windows.
%   db_svVals: Structure of slope variance values, contains basin IDs and 
%       slope variance windows, as well as mean, median, std, min, and max 
%       slope variances. Statistics matrices have rows the same number as 
%       basins, and columns the same number as slope variance windows.

%% Setup
[DBg,~,~] = GRIDobj2mat(DB);
DBg(DBg==0) = NaN;
dbi = unique(DBg(:));
dbi(isnan(dbi)) = [];

db_rVals.BasinIDs = dbi;
db_rVals.Windows = [];
db_rVals.Means = zeros(length(dbi),length(rVals))*NaN;
db_rVals.Medians = db_rVals.Means;
db_rVals.Stds = db_rVals.Means;
db_rVals.Mins = db_rVals.Means;
db_rVals.Maxes = db_rVals.Means;

db_svVals.BasinIDs = dbi;
db_svVals.Windows = [];
db_svVals.Means = zeros(length(dbi),length(svVals))*NaN;
db_svVals.Medians = db_svVals.Means;
db_svVals.Stds = db_svVals.Means;
db_svVals.Mins = db_svVals.Means;
db_svVals.Maxes = db_svVals.Means;

%% Condense roughnesses
rGrids = cell(length(rVals),1);
for i = 1:length(rVals)
    db_rVals.Windows(i) = rVals(i).TrueWindowRes;
    [tmpR,~,~] = GRIDobj2mat(rVals(i).Roughness_Grid);
    rGrids{i} = tmpR;
end

%% Condense slope variance
svGrids = cell(length(svVals),1);
for i = 1:length(svVals)
    db_svVals.Windows(i) = svVals(i).TrueWindowRes;
    [tmpSV,~,~] = GRIDobj2mat(svVals(i).SlopeVariance_Grid);
    svGrids{i} = tmpSV;
end

curPer = 0;
for i = 1:length(dbi)
    if verbose > 0 && i/length(dbi)>=curPer
        disp(sprintf('      %d%% Complete (%d / %d)',round(curPer*100,0),i,length(dbi)))
        curPer = curPer + .1;
    end

    %% Loop through roughness windows
    for j = 1:length(rVals)
        tmpR = rGrids{j};
        tmpR(DBg~=dbi(i)) = NaN;
        tmpR(isnan(tmpR)) = [];

        db_rVals.Means(i,j) = mean(tmpR);
        db_rVals.Medians(i,j) = median(tmpR);
        db_rVals.Stds(i,j) = std(tmpR);
        db_rVals.Mins(i,j) = min(tmpR);
        db_rVals.Maxes(i,j) = max(tmpR);
    end

    %% Loop through slope variance windows
    for j = 1:length(rVals)
        tmpSV = svGrids{j};
        tmpSV(DBg~=dbi(i)) = NaN;
        tmpSV(isnan(tmpSV)) = [];

        if ~isempty(tmpSV)
            db_svVals.Means(i,j) = mean(tmpSV);
            db_svVals.Medians(i,j) = median(tmpSV);
            db_svVals.Stds(i,j) = std(tmpSV);
            db_svVals.Mins(i,j) = min(tmpSV);
            db_svVals.Maxes(i,j) = max(tmpSV);
        else
            db_svVals.Means(i,j) = NaN;
            db_svVals.Medians(i,j) = NaN;
            db_svVals.Stds(i,j) = NaN;
            db_svVals.Mins(i,j) = NaN;
            db_svVals.Maxes(i,j) = NaN;
        end
    end
end
end