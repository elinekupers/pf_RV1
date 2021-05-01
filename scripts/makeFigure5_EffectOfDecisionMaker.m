function [] = makeFigure5_EffectOfDecisionMaker()
% Function to make Figure 5 of the manuscript:
%   Radial asymmetries around the visual field: From retina to cortex to 
%   behavior. By Kupers, Benson, Carrasco, Winawer.
%    YEAR. JOURNAL. DOI.

% This function requires you to have classifier performance on the 
% simulated data. To run this and the RGC model see s_runRGCmodelExample.m


baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
subfolder  = 'onlyL';

plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'idealobserver', subfolder, 'Ideal')
plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', subfolder, 'SVM-Fourier')
plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', subfolder, 'SVM-Energy')
