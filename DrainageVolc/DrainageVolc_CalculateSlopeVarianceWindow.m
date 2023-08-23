function slVarGrid = DrainageVolc_CalculateSlopeVarianceWindow(slp,slpSize)
% Name: DrainageVolc_CalculateSlopeVarianceWindow
% Author: Daniel O'Hara
% Date: 08/18/2023 (mm/dd/yyyy)
% Description: Script to calculate the slope variance (std / mean) of
% pixels over a given window size.
%
% Input:
%   slp: GRIDobj of slope.
%   slpSize: Integer window size.
%
% Output:
%   slVarGrid: GRIDobj of slope variances.

% Get array
[Sg,x,y] = GRIDobj2mat(slp);
[X,Y] = meshgrid(x,y);

% Calculate mean and std arrays
meanS = ones(size(Sg))*NaN;
stdS = meanS;
halfSize = floor(slpSize/2);
for i = 1:size(Sg,1)
    for j = 1:size(Sg,2)
        if isnan(Sg(i,j))
            continue;
        end

        i1 = max([1,i-halfSize]);
        i2 = min([size(Sg,1),i+halfSize]);
        j1 = max([1,j-halfSize]);
        j2 = min([size(Sg,2),j+halfSize]);

        tmp = Sg(i1:i2,j1:j2);
        meanS(i,j) = nanmean(tmp(:));
        stdS(i,j) = nanstd(tmp(:));
    end
end

% Variance
slVar = stdS./meanS;
slVarGrid = GRIDobj(X,Y,slVar);
end