function pack = MorVolc_DefaultVals(pack)
% Name: MorVolc_DefaultVals
% Author: Daniel O'Hara
% Date: 06/30/2021 (mm/dd/yyyy)
% Description: Script to check the fields of the input structure and fill 
%   any missing fields with default values.

%% Check that TopoToolbox is accessible
if ~exist("GRIDobj","file")
    error('This package requires TopoToolbox, downloadable from https://topotoolbox.wordpress.com/.')
end

%% Fill Most Fields if Missing
packNumFields = {'contIter';...
    'dx';'hypsIter';'peakDiff';'peakContIter';...
    'craterContIter';'summitRegion';...
    'correctCrater';'ignoreSummitIndex';...
    'plotResults';'verbose';'saveInputs';...
    'zipFiles';'deleteAfterZip';'visPlots';'xlsFile'};

packEmptyFields = {'maskXY';'craterXY';...
    'saveFigFolder';'saveResFolder';'figTitlePrefix'};

packStringFields = {'figPrefix'};

folderPathFields = {'saveFigFolder';'saveResFolder'};

packNumVals = {-.1;...
    30;-.01;0.05;-.1;...
    -.1;1;...
    1;1;...
    0;1;0;...
    0;0;1;0};

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
        eval(sprintf('pack.%s = '''';',packStringFields{i}))
    end
end

%% Fill Interpolation Fields if Missing
if ~isfield(pack,'interpSurfaces')
    disp('WARNING: Missing field ''interpSurfaces'', filling with default value.')
    pack.interpSurfaces.Natural = 1;
    pack.interpSurfaces.IDW = 0;
    pack.interpSurfaces.Kriging = 0;
else
    if ~isfield(pack.interpSurfaces,'Natural')
        disp('WARNING: Missing field ''interpSurfaces.Natural'', filling with default value.')
        pack.interpSurfaces.Natural = 1;
    end
    
    if ~isfield(pack.interpSurfaces,'IDW')
        disp('WARNING: Missing field ''interpSurfaces.IDW'', filling with default value.')
        pack.interpSurfaces.IDW = 0;
    end
    
    if ~isfield(pack.interpSurfaces,'Kriging')
        disp('WARNING: Missing field ''interpSurfaces.Kriging'', filling with default value.')
        pack.interpSurfaces.Kriging = 0;
    end
end

%% Fill Roughness Windows if Missing
if ~isfield(pack,'roughnessWindows')
    disp('WARNING: Missing field ''roughnessWindows'', filling with default value.')
    pack.roughnessWindows = [250,500,1000,2000];
end

if ~isfield(pack,'roughnessType')
    disp('WARNING: Missing field ''roughnessType'', filling with default value.')
    pack.roughnessType = 'tpi';
end

if ~isfield(pack,'slopeVarianceWindows')
    disp('WARNING: Missing field ''slopeVarianceWindows'', filling with default value.')
    pack.slopeVarianceWindows = [250,500,1000,2000];
end

%% Test folder paths
for i = 1:length(folderPathFields)
    if ~isempty(folderPathFields{i})
        evalc(sprintf('pack.%s = strrep(pack.%s,''\\'',''/'');',folderPathFields{i},folderPathFields{i}));
        evalc(sprintf('endChar = pack.%s(end);',folderPathFields{i}));
        if ~strcmp(endChar,'/')
            evalc(sprintf('pack.%s(end) = ''/''',folderPathFields{i}));
        end
    end
end

%% Fill version
pack.version = '9.08192025';