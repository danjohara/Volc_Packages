function kn_ID_BID_XY_dist_dz = DrainageVolc_Find_Knickpoints(DEM,DB,channelThreshold,tol,parallelProc,verbose)
% Name: DrainageVolc_Find_Knickpoints
% Author: Daniel O'Hara
% Date: 06/09/2024 (mm/dd/yyyy)
% Description: Script to determine knickpoints on the edifice.
%
% Input:
%   DEM: Cut GRIDobj of DEM.
%   DB: GRIDobj of drainage basins.
%   channelThreshold: Drainage area threshold for channelization.
%   tol: Channel profile tolerance for knickpoints.
%   parallelProc: Flag for whether to use parallelization.
%   verbose: Flag for whether updates should be outputted.
%   
% Output:
%   kn_ID_BID_XY_dist_dz: Array of knickpoint values. Contains a knickpoint
%       unique ID, associated basin ID, knickpoint X-Y coordinates,
%       knickpoint upstream distance, and the magnitude of the knickpoint
%       give as the elevation difference.

kn_ID_BID_XY_dist_dz = [];
counter = 1;

uniDB = unique(DB.Z(:));
uniDB(isnan(uniDB)) = [];

warning('off','all')

curPer = 0;
for i = 1:length(uniDB)
    if verbose > 0 && i/length(uniDB)>=curPer
        disp(sprintf('         %d%% Complete (%d / %d)',round(curPer*100,0),i,length(uniDB)))
        curPer = curPer + .1;
    end

    %% Isolate Basin
    tmpDEM = DEM;
    tmpDEM.Z(DB.Z~=uniDB(i)) = NaN;

    if sum(~isnan(tmpDEM.Z(:)))*tmpDEM.cellsize^2 < channelThreshold
        continue;
    end
    
    tmpDEM = crop(tmpDEM);

    %% Get Basin STREAMobj
    FD = FLOWobj(tmpDEM);
    S = STREAMobj(FD,'minarea',channelThreshold,'unit','mapunits');

    if isempty(S.ix)
        continue;
    end

    %% Find Knickpoints
    [~,kp] = knickpointfinder(S,tmpDEM,'split',parallelProc,'tol',tol,'plot',false,'verbose',false);

    if kp.n == 0
        continue;
    end

    for j = 1:kp.n
        kn_ID_BID_XY_dist_dz = [kn_ID_BID_XY_dist_dz;counter,uniDB(i),kp.x(j),kp.y(j),kp.distance(j),kp.dz(j)];
        counter = counter + 1;
    end
end
warning('on','all')

end