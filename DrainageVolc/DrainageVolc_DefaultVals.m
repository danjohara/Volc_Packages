function pack = DrainageVolc_DefaultVals(pack)
% Name: Drainagevolc_DefaultVals
% Author: Daniel O'Hara
% Date: 06/30/2021 (mm/dd/yyyy)
% Description: Script to check the fields of the input structure and fill 
%   any missing fields with default values.

%% Check that TopoToolbox is accessible
if ~exist("GRIDobj","file")
    error('This package requires TopoToolbox, downloadable from https://topotoolbox.wordpress.com/.')
end

%% Fill Most Fields if Missing
packNumFields = {'dx';'hypsIter';'channelThreshold';...
    'basinTopN';'basinTopNFromContLength';'parallelProc';...
    'plotResults';'dynamicThreshold';'dynamicThresholdPixelStep';...
    'Analyze_Divides';...
    'basinContIter';'smoothBasinPointWavelength';...
    'MN';'limitHacksLaw';'slopeVariance_ContIter';'basinRadIter'};

packEmptyFields = {'craterXY';'maskMap';'saveFigFolder';'saveResFolder'};

packStringFields = {'figPrefix';'roughnessType'};

packNumVals = {30;.01;5e5;...
    -.3;0;0;...
    0;1;20;...
    0;...
    -.05;1000;...
    NaN;0;-.05;-.01};

packStringVals = {'''';'tpi'};

for i = 1:length(packNumFields)
    if ~isfield(pack,packNumFields{i})
        disp(sprintf('WARNING: Missing field ''%s'', filling with default value.',packNumFields{i}))
        eval(sprintf('pack.%s = %f;',packNumFields{i},packNumVals{i}))
    end
end

for i = 1:length(packEmptyFields)
    if ~isfield(pack,packEmptyFields{i})
        disp(sprintf('WARNING: Missing field ''%s'', filling with default value.',packEmptyFields{i}))
        eval(sprintf('pack.%s = [];',packEmptyFields{i}))
    end
end

for i = 1:length(packStringFields)
    if ~isfield(pack,packStringFields{i})
        disp(sprintf('WARNING: Missing field ''%s'', filling with default value.',packStringFields{i}))
        eval(sprintf('pack.%s = ''%s'';',packStringFields{i},packStringVals{i}))
    end
end

%% Fill Roughness Windows if Missing
if ~isfield(pack,'roughnessWindows')
    disp('WARNING: Missing field ''roughnessWindows'', filling with default value.')
    pack.roughnessWindows = [250,500,1000,2000];
end

if ~isfield(pack,'slopeVarianceWindows')
    disp('WARNING: Missing field ''slopeVarianceWindows'', filling with default value.')
    pack.slopeVarianceWindows = [250,500,1000,2000];
end

%% Version Number
pack.version = '6.08202023';