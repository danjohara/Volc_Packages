% Name: TestMorVolc
% Author: Daniel O'Hara
% Date: 02/23/2021 (mm/dd/yyyy)
% Description: Example script to test the MorVolc_Analysis function
%   for the Aracar volcano.

%% Script Parameters
packagePath = '.\..';
tifFolder = '.\..\Example_DEMs\Aracar_Data\';
tifFile = 'Aracar_Elevation_UTM.tif';
boundaryFile = 'Aracar_Boundary.shp';
craterFile = 'Aracar_Crater.shp';
maskFile = 'Aracar_Mask.shp';
dx = 30;
prefix = 'Aracar_';

%% Create Package
pack.tifFile = [tifFolder,tifFile];
pack.boundaryXY = [tifFolder,boundaryFile];
pack.maskMap = []; % [tifFolder,maskFile]
pack.craterXY = [tifFolder,craterFile]; % [tifFolder,craterFile]

pack.dx = dx;
pack.peakDiff = .05;
pack.contIter = -.05;
pack.peakContIter = -.05;
pack.craterContIter = -.05;

pack.correctCrater = 1;
pack.summitRegion = 1;
pack.ignoreSummitIndex = 1;

pack.interpSurfaces.Natural = 1;
pack.interpSurfaces.IDW = 1;
pack.interpSurfaces.Kriging = 1;

pack.xlsFile = '.\..\Example_Results\Morvolc\Test_Results.xlsx';
pack.saveResFolder = '.\..\Example_Results\MorVolc\';

pack.plotResults = 1;
pack.figPrefix = prefix;
pack.saveFigFolder = '.\..\Example_Results\MorVolc\';

%% Run Analysis
addpath(genpath(packagePath))
morRes = MorVolc_Analysis(pack);