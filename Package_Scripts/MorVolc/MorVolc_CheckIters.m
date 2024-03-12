function MorVolc_CheckIters(contIter,craterXYZ,craterContIter,peakContIter)
% Name: MorVolc_CheckIters
% Author: Daniel O'Hara
% Date: 06/29/2022 (mm/dd/yyyy)
% Description: Script to input values are correct.

if contIter == 0 || contIter <= -1
    error('Invalid contIter; value needs to be [-1,0] or > 0');
end

if ~isempty(craterXYZ) && (craterContIter == 0 || craterContIter <= -1)
    error('Invalid craterContIter; value needs to be [-1,0] or > 0');
end

if peakContIter == 0 || peakContIter <= -1
    error('Invalid peakContIter; value needs to be [-1,0] or > 0');
end
