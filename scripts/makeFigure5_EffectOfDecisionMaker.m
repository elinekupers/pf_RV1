function [] = makeFigure5_EffectOfDecisionMaker()

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
subfolder  = 'onlyL';

plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'idealobserver', subfolder, 'Ideal')
plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', subfolder, 'SNR')
plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', subfolder, 'SVM')