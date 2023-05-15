function [basinAreaThreshold,tranStartA,MN] = DrainageVolc_SA_PiecewiseRegression(topNFlowAS,areaStep)
% Name: SA_PiecewiseRegression
% Author: Daniel O'Hara
% Date: 07/09/2021 (mm/dd/yyyy)
% Description: Script to find the best drainage area threshold for
%   channelization using piece-wise regression to determine drainage area
%   that gives breaks in slope-area relationship.
%
% Input:
%   TopNFlowAS: Drainage area and slope values of the flow paths associated
%       with the largest basins.
%   areaStep: Drainage area step value to find break in relationship.
%
% Output:
%   basinAreaThreshold: Array containing best-fit drainage area channel
%       threshold, r^2 value of the first (hillslope) regresssion, r^2
%       value of the second (channel) regression, and the r' value given by
%       r' = (|r_1 - r_2|) / (r_1 + r_2).
%   tranStartA: Area threhold for transition zone 
%       initiation.

posSTest = min(topNFlowAS(:,1))+areaStep:areaStep:max(topNFlowAS(:,1))-areaStep;
posS = zeros(size(posSTest));

for i = 1:size(posSTest,1)
    AS = topNFlowAS(topNFlowAS(:,1)<=posSTest(i),:);

    try
        mdl = fitlm(log10(AS(:,1)),log10(AS(:,2)));
    catch er
        posS(i) = NaN;
        continue;
    end
    
    posS(i) = mdl.Coefficients.Estimate(2);
end

ii = find(posS>0,1,'last');
if isempty(ii)
    tranStartA = posSTest(1);
else
    tranStartA = posSTest(ii);
end
        

% areaRTests = min(topNFlowAS(:,1))+areaStep:areaStep:max(topNFlowAS(:,1))-areaStep;
areaRTests = tranStartA+areaStep:areaStep:max(topNFlowAS(:,1))-areaStep;
areaRTests = [areaRTests',ones(length(areaRTests),3)];
topNFlowAS(topNFlowAS(:,2)==0,:) = [];
topNFlowAS(:,2) = tand(topNFlowAS(:,2));

for i = 1:size(areaRTests,1)
    AS1 = topNFlowAS(topNFlowAS(:,1)<=areaRTests(i,1),:);
    AS2 = topNFlowAS(topNFlowAS(:,1)>areaRTests(i,1),:);
    
    try
        mdl1 = fitlm(log10(AS1(:,1)),log10(AS1(:,2)));
        mdl2 = fitlm(log10(AS2(:,1)),log10(AS2(:,2)));
    catch er
        areaRTests(i,2:4) = [1,1,1]*NaN;
        continue;
    end
    
    r2_1 = mdl1.Rsquared.Ordinary;
    r2_2 = mdl2.Rsquared.Ordinary;
    
    areaRTests(i,2:3) = [r2_1,r2_2];
%     areaRTests(i,4) = (r2_1+r2_2)/2;
%     areaRTests(i,4) = ((r2_1+r2_2)+abs(r2_1-r2_2))/(r2_1+r2_2);
%     areaRTests(i,4) = (abs(r2_1-r2_2))/(r2_1+r2_2); % previous chosen one.
    areaRTests(i,4) = min([r2_1,r2_2])/max([r2_1,r2_2]);
end
areaRTests(areaRTests(:,4)==0,2:4) = NaN;
% Get best minA based on balance between R1 & R2.
areaRTests(areaRTests(:,1)>5e6,:) = NaN;
minA_i = find(areaRTests(:,4)==nanmax(areaRTests(:,4)),1); % use nanmin if previous

% get best minA based on best R2 within a 5 km^2 limit.
% areaRTests(areaRTests(:,1)>5e6,:) = NaN;
% minA_i = find(areaRTests(:,3)==nanmax(areaRTests(:,3)),1);

% Continue rest of script.
basinAreaThreshold = areaRTests(minA_i,:);

AS = topNFlowAS(topNFlowAS(:,1)>basinAreaThreshold(1),:);
mdl = fitlm(log10(AS(:,1)),log10(AS(:,2)));
MN = mdl.Coefficients.Estimate(2);
end