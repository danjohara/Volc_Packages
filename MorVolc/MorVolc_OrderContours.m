function allCC_XY = MorVolc_OrderContours(ccXY)
%%
% Name: MorVolc_OrderContours
% Date: 08/23/2023 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to spatially order contour lines
%
% Input:
%   ccXY: Cell array of contour lines.
%
% Output:
%   allCC_XY: 2xN array of ordered contour XY data.

sXY = zeros(length(ccXY),2);
fXY = sXY;
for i = 1:length(ccXY)
    sXY(i,:) = ccXY{i}(1,:);
    fXY(i,:) = ccXY{i}(end,:);
end

sfXY = [sXY;fXY];
sfDists = triu(pdist2(sfXY,sfXY));

disp('here')
[ii,jj] = find(sfDists == max(sfDists(:)),1);

allCC_XY = [];
tmpSXY = sXY;
tmpFXY = fXY;
tmpCCXY = ccXY;
if ii <= size(sXY,1)
    allCC_XY = ccXY{ii};
    cutI = ii;
else
    allCC_XY = flipud(ccXY{ii-size(sXY,1)});
    cutI = ii-length(sXY);
end

tmpSXY(cutI,:) = [];
tmpFXY(cutI,:) = [];
tmpCCXY(cutI) = [];

while ~isempty(tmpFXY)
    tmpSFXY = [tmpSXY;tmpFXY];
    sfDists = pdist2(allCC_XY(end,:),tmpSFXY);
    
    ii = find(sfDists==min(sfDists),1);

    if ii <= size(tmpSXY,1)
        allCC_XY = [allCC_XY;NaN,NaN;tmpCCXY{ii}];
        cutI = ii;
    else
        allCC_XY = [allCC_XY;NaN,NaN;flipud(tmpCCXY{ii-size(tmpSXY,1)})];
        cutI = ii-size(tmpSXY,1);
    end

    tmpSXY(cutI,:) = [];
    tmpFXY(cutI,:) = [];
    tmpCCXY(cutI) = [];
end
end