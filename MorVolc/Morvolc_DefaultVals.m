function pack = Morvolc_DefaultVals(pack)
% Name: Morvolc_DefaultVals
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
    'dx';'peakDiff';'peakContIter';...
    'craterContIter';'summitRegion';...
    'correctCrater';'ignoreSummitIndex';...
    'plotResults'};

packEmptyFields = {'maskMap';'craterXY';...
    'saveFigFolder';'xlsFile';'saveResFolder'};

packStringFields = {'figPrefix'};

packNumVals = {-.1;...
    30;0.05;-.1;...
    -.1;1;...
    1;1;...
    0};

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

%% Fill version
pack.version = '6.08202023';