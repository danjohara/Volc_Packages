function db_contSin = DrainageVolc_Collect_DB_ContourSinuosity(DEM,DB,cIter,verbose)
% Name: DrainageVolc_Collect_DB_ContourSinuosity
% Author: Daniel O'Hara
% Date: 06/07/2024 (mm/dd/yyyy)
% Description: Script to calculate sinuosity values of contours along each
%   basin, as a representation of incision.
%
% Input:
%   DEM: GRIDobj of elevations.
%   DB: GRIDobj of drainage basins.
%   cIter: Contour interval, can be given as a whole number or as a
%       percentage (with a negative number between 0 and -1) of the entire
%       range of elevations.
%
% Output:
%   db_contSin: Structure of basin contour sinuosities. Contains analyzed
%       basin IDs, contour elevatoins, normalized contour elevations, 
%       array of contour sinuosities (rows are basins, columns are 
%       contours), as well as mean, median, std, min, and max sinuosity of 
%       each basin. 

[Zg,x,y] = GRIDobj2mat(DEM);
x = double(x);
y = double(y);
Zg = double(Zg);

[DBg,~,~] = GRIDobj2mat(DB);
DBg(isnan(Zg)) = NaN;
dbi = unique(DBg(:));
dbi(isnan(dbi)) = [];

if cIter > -1 && cIter < 0
    cIter = abs(cIter)*range(Zg(:));
end

conts = min(Zg(:))+cIter:cIter:max(Zg(:));
normConts = conts-min(Zg(:));
normConts = normConts./max(normConts);

db_contSin.BasinIDs = dbi;
db_contSin.Contours = conts;
db_contSin.Norm_Contours = normConts;
db_contSin.Individual_Sinuosities = zeros(length(dbi),length(conts))*NaN;

curPer = 0;
for i = 1:length(dbi)
    if verbose > 0 && i/length(dbi)>=curPer
        disp(sprintf('      %d%% Complete (%d / %d)',round(curPer*100,0),i,length(dbi)))
        curPer = curPer + .1;
    end

    %% Isolate topography
    tmpZ = Zg;
    tmpZ(DBg~=dbi(i)) = NaN;

    %% Loop through contours
    for j = 1:length(conts)
        cc = contourc(x,y,tmpZ,[1,1]*conts(j));

        if isempty(cc)
            continue;
        end

        cc = Convert_Contours(cc,1,0);
        contPerm = sum(sqrt(sum(diff(cc,1,1).^2,2)));
        eucDist = sqrt((cc(1,1)-cc(end,1))^2 + (cc(1,2)-cc(end,2))^2);

        db_contSin.Individual_Sinuosities(i,j) = contPerm/eucDist;
    end
end

db_contSin.Means = nanmean(db_contSin.Individual_Sinuosities,2);
db_contSin.Medians = nanmedian(db_contSin.Individual_Sinuosities,2);
db_contSin.Stds = nanstd(db_contSin.Individual_Sinuosities,0,2);
db_contSin.Mins = nanmin(db_contSin.Individual_Sinuosities,[],2);
db_contSin.Maxes = nanmax(db_contSin.Individual_Sinuosities,[],2);
end