function [tmpPF0,PF,lDev_Map,lDev_Titles] = DrainageVolc_HacksLaw(DB,BID,BA,BL,channelThreshold)
% Name: DrainageVolc_HacksLaw
% Author: Daniel O'Hara
% Date: 07/03/2023 (mm/dd/yyyy)
% Description: Script to calculate the power-law function between basin
%   areas and lengths (i.e., Hack's Law).
%
% Input:
%   DB: GRIDobj of drainage basins.
%   BID: Nx1 array of basin IDs.
%   BA: Nx1 array of basin areas.
%   BL: Nx1 array of basin lengths.
%   channelThreshold: Drainage area threshold for channelization (NaN if 
%       not using this constraint).
% Output:
%   tmpPF0: Matrix of basin IDs, basin areas, basin lengths, basin length
%       deviation from the established power-law, and flag for whether the
%       basin was used to derive the power law.
%   PF: Power-law parameters (PF(1)*(basinArea)^(PF(2)).
%   lDev_Map: GRIDobj of basin length deviation from the power law.
%   lDev_Titles: Titles associated with tmpPF0.

%% Get Basin Data
tmpPF = [BID,BA,BL];
tmpPF0 = [tmpPF,zeros(size(tmpPF(:,1))),ones(size(tmpPF(:,1)))];
if ~isnan(channelThreshold)
    tmpPF0(tmpPF(:,2)<channelThreshold,5) = 0;
    tmpPF(tmpPF(:,2)<channelThreshold,:) = [];
end

%% Calculate Power-Law Fit
tmpPF(sum(isnan(tmpPF),2)>0,:) = [];
tmpPF(tmpPF(:,3)==0,:) = [];
powerFit = polyfit(log10(tmpPF(:,2)),log10(tmpPF(:,3)),1);
PF = [10^powerFit(2),powerFit(1)];

%% Deviation from Power-Law
tmpPF0(:,4) = log10(tmpPF0(:,3)) - log10(PF(1)*tmpPF0(:,2).^PF(2));
lDev_Map = DB;
lDev_Map.Z(:) = NaN;
for i = 1:size(tmpPF0,1)
    if tmpPF0(i,5)==1
        lDev_Map.Z(DB.Z==tmpPF0(i,1)) = tmpPF0(i,4);
    end
end

lDev_Titles = {'Basin ID','Drainage Area','Length','Deviation','Used in Power Fit'};

end