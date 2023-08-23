function [DBxy,DBxy_cell] = DrainageVolc_GetBasinBoundaries(DB,A,dx,channelThreshold)
% Name: GetBasinBoundaries
% Author: Daniel O'Hara
% Date: 06/15/2021 (mm/dd/yyyy)
% Description: Script to collect the boundaries of all relevant basins.
%   Relevant basins are considered to be those have minimum drainage area
%   equivelant to the channel threshold.
%
% Input:
%   DB: GRIDobj of drainage basins.
%   A: GRIDobj of drainage areas (assumed to be in number of pixels).
%   dx: Grid resolution.
%   channelThreshold
%
% Output:
%   DBxy: N x 2 array of drainage basin X- and Y- coordinates, separated by
%       NaN's for each basin.
%   DBxy_cell: Cell array of drainage basin X- and Y- coordinates.
DBxy = [];
DBxy_cell = {};
[DBg,X,Y] = GRIDobj2mat(DB);
[Ag,~,~] = GRIDobj2mat(A);
Ag = Ag*dx^2;
dbu = unique(DBg(:));
dbu(isnan(dbu)) = [];

for i = 1:length(dbu)
    tmpAg = Ag;
    tmpAg(DBg~=dbu(i)) = 0;
    if nanmax(tmpAg(:)) < channelThreshold
        continue
    end
    
    bb = bwboundaries(DBg==dbu(i));
    tmp = [];
    for j= 1:size(bb{1},1)
        tmp = [tmp;X(bb{1}(j,2)),Y(bb{1}(j,1))];
    end
    DBxy = [DBxy;tmp;NaN,NaN];
    DBxy_cell = [DBxy_cell;{tmp}];
end