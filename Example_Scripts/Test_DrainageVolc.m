% Name: Test_DrainageVolc
% Author: Daniel O'Hara
% Date: 02/25/2021 (mm/dd/yyyy)
% Description: Example script to run the DrainageVolc algorithm
%   for the Aracar volcano.

%% Script Parameters
packagePath = './..';
tifFolder = './../Example_DEMs/Aracar_Data/';
tifFile = 'Aracar_Elevation_UTM.tif';
boundaryFile = 'Aracar_Boundary.shp';
craterFile = 'Aracar_Crater.shp';
maskFile = 'Aracar_Mask.shp';
dx = 30;
prefix = 'Aracar_Test_';

%% Create Package
pack.tifFile = [tifFolder,tifFile];
pack.boundaryXY = [tifFolder,boundaryFile];
pack.maskMap = []; 
% pack.maskMap = [tifFolder,maskFile] % Uncomment to mask sections of the volcano.
pack.craterXY = [tifFolder,craterFile];

pack.dx = dx;
pack.hypsIter = .01;
pack.roughnessWindows = [250,500,1000,2000];

pack.channelThreshold = 5e5;
pack.dynamicThreshold = 1;
pack.dynamicThresholdPixelStep = 20;

pack.basinContIter = -.05;
pack.basinTopN = -.3; 
pack.basinTopNFromContLength = 1;

pack.smoothBasinPointWavelength = 1000;
pack.parallelProc = 0;
pack.Analyze_Divides = 1;
pack.MN = NaN;

pack.saveResFolder = './../Example_Results/DrainageVolc/';

pack.plotResults = 1;
pack.figPrefix = prefix;
pack.saveFigFolder = './../Example_Results/DrainageVolc/';

%% Run Analysis
addpath(genpath(packagePath))
mets = DrainageVolc_Analysis(pack);