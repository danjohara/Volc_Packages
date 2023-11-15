function C = Convert_Contours(cc,onlyLargest,withElevations)
%%
% Name: Convert_Contours
% Author: Daniel O'Hara
% Data: 03/17/2021 (mm/dd/yyyy)
% Description: Function to convert the matrix format of contours from the
%   contourc function into a cell array of individual contours. If
%   onlyLargest is set to 1, the script will only return the largest
%   contour, otherwise it will return all.
%
% Input:
%   cc: Matrix format of contours, as given by contourc/contourf/etc.
%   onlyLargest: Flag to return only the largest contour
%
% Output: 
%   C: Cell array of individual contour coordinates, or array of largest 
%       contour coordinates.

C = {};

maxIter = NaN;
maxSize = -Inf;

i = 1;
while i<size(cc,2)
    if nargin == 3 && withElevations == 1
        tmp = [cc(:,i+1:i+cc(2,i))]';
        tmp = [tmp,ones(size(tmp(:,1)))*cc(1,i)];
        C = [C;{tmp}];
    else
        C = [C;{[cc(:,i+1:i+cc(2,i))]'}];
    end
    i = i+cc(2,i)+1;

    if size(C{end},1) > maxSize
        maxSize = size(C{end,1});
        maxIter = length(C);
    end
end

if onlyLargest && ~isnan(maxIter)
    C = C{maxIter};
end
end