% Name: TestMorVolc
% Author: Daniel O'Hara
% Date: 02/23/2021 (mm/dd/yyyy)
% Description: Example script to test the MorVolc_Analysis function
%   for the Aracar volcano.

%% Script Parameters
packagePath = './../..';
tifFolder = './../Example_DEMs/Aracar_Data/';
tifFile = 'Aracar_Elevation_UTM.tif';
boundaryFile = 'Aracar_Boundary.shp';
craterFile = 'Aracar_Crater.shp';
maskFile = 'Aracar_Mask.shp';
dx = 30;
prefix = 'Aracar_';

%% Create Package
% Input Files
pack.tifFile = [tifFolder,tifFile];
pack.boundaryXY = [tifFolder,boundaryFile];
pack.maskXY = []; % [tifFolder,maskFile]
pack.craterXY = [tifFolder,craterFile]; % [tifFolder,craterFile]

% Map Resolution
pack.dx = dx;

% Analysis elevation intervals
pack.peakDiff = .05;
pack.contIter = -.05;
pack.peakContIter = -.05;
pack.craterContIter = -.05;

% Summit analysis
pack.correctCrater = 1;
pack.summitRegion = 1;
pack.ignoreSummitIndex = 1;

% Basal surface interpolation
pack.interpSurfaces.Natural = 1;
pack.interpSurfaces.IDW = 1;
pack.interpSurfaces.Kriging = 1;

% Topography analysis
pack.hypsIter = .01;
pack.roughnessWindows = [250,500,1000,2000];
pack.roughnessType = 'tpi';
pack.slopeVarianceWindows = [250,500,1000,2000];

% Save parameters
pack.xlsFile = 1;
pack.saveResFolder = './../Example_Results/MorVolc/';
pack.saveFigFolder = './../Example_Results/MorVolc/';
pack.saveInputs = 1;

% Figure parameters
pack.plotResults = 1;
pack.visPlots = 1;
pack.figTitlePrefix = 'Aracar';
pack.figPrefix = prefix;

% Verbose
pack.verbose = 1;

% Zipping and deleting results
pack.zipFiles = 3;
pack.deleteAfterZip = 0;

%% Run Analysis
addpath(genpath(packagePath))
morRes = MorVolc_Analysis(pack);