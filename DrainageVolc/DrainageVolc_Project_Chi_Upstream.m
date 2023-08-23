function Upstream_Chi = DrainageVolc_Project_Chi_Upstream(DEM,DB,chiS,Chi)
%%
% Name: DrainageVolc_Project_Chi_Upstream
% Date: 09/01/2022
% Author: Daniel O'Hara
% Description: Script to project channel \Chi values upstream in
%   topographic order to the divides.
%
% Input: 
%   DEM: GRIDobj of elevations.
%   DB: GRIDobj of basins.
%   chiS: STREAMobj of channels with \Chi values.
%   Chi: Array of \Chi values.
%
% Output:
%   Upstream_Chi: GRIDobj of upstream-projected \Chi values.

%% Get Grids
[DBg,~,~] = GRIDobj2mat(DB);
[Zg,~,~] = GRIDobj2mat(DEM);

%% Collect all basin IDs, elevations, and chi values of stream network
IX_dbID_Z_Chi = zeros(size(Chi,1),4)*NaN;

for i = 1:length(chiS.IXgrid)
    IX_dbID_Z_Chi(i,:) = [chiS.IXgrid(i),DBg(chiS.IXgrid(i)),Zg(chiS.IXgrid(i)),Chi(i)];
end

%% Create upstream chi grid and fill values
Upstream_Chi = DEM;
Upstream_Chi.Z(:) = NaN;

db_uni = unique(IX_dbID_Z_Chi(:,2));
for i = 1:length(db_uni)
    basinTest = DBg==db_uni(i);
    
    tmp_IX_dbID_Z_Chi = IX_dbID_Z_Chi;
    tmp_IX_dbID_Z_Chi(tmp_IX_dbID_Z_Chi(:,2)~=db_uni(i),:) = [];
    tmp_IX_dbID_Z_Chi = sortrows(tmp_IX_dbID_Z_Chi,3,"ascend");

    for j = 1:size(tmp_IX_dbID_Z_Chi,1)
        Ztest = Zg>=tmp_IX_dbID_Z_Chi(j,3);
        Upstream_Chi.Z((basinTest.*Ztest) == 1) = tmp_IX_dbID_Z_Chi(j,4);
    end
end