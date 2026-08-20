% Name: Test_DrainageVolc
% Author: Daniel O'Hara
% Date: 02/25/2021 (mm/dd/yyyy)
% Description: Example script to run the DrainageVolc_Analysis function
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
pack.craterXY = [tifFolder,craterFile]; %[tifFolder,craterFile];

% Map Resolution
pack.dx = dx;

% Topography analysis
pack.hypsIter = .01;
pack.roughnessWindows = [250,500,1000,2000];
pack.roughnessType = 'tpi';
pack.slopeVarianceWindows = [250,500,1000,2000];

% Basin analysis
pack.basinContIter = -.05;
pack.basinRadIter = -.01;
pack.basinTopN = -.3; 
pack.basinTopNFromContLength = 1;
pack.limitHacksLaw = 1;
pack.smoothBasinPointWavelength = 1000;
pack.parallelProc = 0;
pack.basinStatThreshold = 1e5;
pack.contourSinuosity_ContIter = -.1;
pack.fillDEM = 1;

% Channel values
pack.channelThreshold = 5e5;
pack.dynamicThreshold = 0;
pack.dynamicThresholdPixelStep = 20;
pack.conformityWavelength = NaN;
pack.conformityStreamDist = 500;
pack.MN = NaN;
pack.knickpointTolerance = 5;
pack.chi_Zcutoff = NaN;
pack.chi_removeUpperBasins = 0;
pack.concavityType = 'lad';

% Divide Analysis
pack.Analyze_Divides = 0;
pack.Divide_Order_Cutoff = 0;
pack.DAI_Integral_BinWidth = .1;

% Saving Results
pack.saveResFolder = './../Example_Results/DrainageVolc/';
pack.saveFigFolder = './../Example_Results/DrainageVolc/';
pack.saveInputs = 1;

% Plotting
pack.plotResults = 1;
pack.visPlots = 1;
pack.figTitlePrefix = 'Aracar';
pack.figPrefix = prefix;
pack.visPlots = 1;

% Verbose:
pack.verbose = 1;

% Zipping and deleting results
pack.zipFiles = 3;
pack.deleteAfterZip = 0;

%% Run Analysis
addpath(genpath(packagePath))
mets = DrainageVolc_Analysis(pack);