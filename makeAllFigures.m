%% makeAllFigures
%
% Make all manuscript figures from the paper:
% Asymmetries around the visual field: From retina to cortex to behavior.
% By Kupers, Benson, Carrasco, Winawer.
% JOURNAL. YEAR. DOI.

% Figure 1: Cone density, mRGC density, V1 CMF vs eccentricity, averaged across meridians
makeFigure1_DensityVSEccentricity

% Figure 2: Cone density, mRGC density, V1 CMF vs eccentricity, as a function of polar angle
makeFigure2_DensityVSPolarAngle

% Figure 3: Overview over computational observer model
makeFigure3_ModelOverview

% Figure 4: RGC layer properties
makeFigure4_RGCRFs

% Figure 5: Effect of  3 different inference engines on RGC responses
makeFigure5_EffectOfDecisionMaker

% Figure 6: For each simulated mRGC:cone ratio show thresholds vs cone 
% density + all at once in 3d mesh
stimTemplateFlag = false; saveFigs = false;
makeFigure6_3DThresholdDensityRatio(stimTemplateFlag,saveFigs);

% Figure 7: Contrast sensitivity, and polar angle asymmetries for cones,
% RGCs and behavior.
makeFigure7_Sensitivity_and_Asymmetry

%% SUPPLEMENTARY FIGURS
% Supplementary Figure 1: Polar angle asymmetries for the different
% computational toolboxes and datsets
makeSupplementalFig1_HVAVMAvsEccen

% Supplementary Figure 2: Can be recomputed by makeFigure1_DensityVSEccentricity

% Supplementary Figure 3: 
% Psychometric functions for RGC responses, ratio 1-5
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
for ratio = 1:5
    plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',ratio)
end

% Supplementary Figure 4: 
% Psychometric functions for cone absorptions + thresholds
plotPsychometricFunctions('conedensity', 'subFolderName', 'average', 'saveFig', false, 'plotAvg', true, 'stimTemplateFlag', false)

% Supplementary Figure 5: 
% Psychometric functions for RGC responses using SVM-Energy template observer
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
for ratio = 1:5
    plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',1, 'stimTemplateFlag', true)
end

% Supplementary Figure 6:
% Thresholds vs cone density, for each simulated mRGC:cone ratio + 3d mesh
stimTemplateFlag = true; saveFigs = false;
makeFigure6_3DThresholdDensityRatio(stimTemplateFlag,saveFigs)

% Supplementary Figure 7:
% Psychometric functions for cone absorptions + thresholds classifier with
% SVM-Energy template decision maker
plotPsychometricFunctions('conedensity', 'subFolderName', 'average_svmEnergy', 'saveFig', false, 'plotAvg', true, 'stimTemplateFlag', true)

% Supplementary Figure 8:
% Psychometric functions for cone current + thresholds classifier with
% SVM-Fourier decision maker
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
plotPsychometricFunctionsConeCurrent(baseFolder, 'conedensity','subFolder','average','plotAvg',true, 'fitTypeName', 'lowess-mesh')
