function DrainageVolc_CheckIters(hypsIter,basinContIter,basinRadIter,basinTopN)
% Name: DrainageVolc_CheckIters
% Author: Daniel O'Hara
% Date: 08/21/2023 (mm/dd/yyyy)
% Description: Script to input values are correct.

if hypsIter <= 0 || hypsIter >= 1
    error('Invalid hypsIter; value needs to be [0,1]');
end

if basinContIter == 0 || basinContIter <= -1
    error('Invalid basinContIter; value needs to be [-1,0] or > 0');
end

if basinRadIter == 0 || basinRadIter <= -1
    error('Invalid basinRadIter; value needs to be [-1,0] or > 0');
end

if basinTopN <= -1
    error('Invalid basinTopN; value needs to be > -1');
end
end