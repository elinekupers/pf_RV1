%% s_runRGCmodelExamples

% A guide to how you can run different types of observer model simulations

%% 1. Stimulus --> Cone absorptions --> Classify with SVM-Fourier

% 1.1 Vary cone density, run 1...5
for runNr = 1:5
    saveFolder = sprintf('run%d',runNr);
    runComputationalObserverModel('conedensity', 'saveFolder', saveFolder, 'seed', runNr)
end

% 1.2 Simple scenario, L-cone only, no eye movements nor stimulus phase
% randomization, run 1
runComputationalObserverModel('defaultnophaseshift', 'saveFolder', 'onlyL', 'seed', 1)

% 1.3 Simple scenario, L-cone only, no eye movements nor stimulus phase
% randomization, run 1, ideal observer
runComputationalObserverModel('idealobserver', 'saveFolder', 'onlyL', 'seed', 1)


% 1.3 Get cone responses for energy template
runComputationalObserverModel('template','saveFolder','run1','seed',1)

% 1.4 Get cone current responses, vary cone density.
for runNr = 1:5
    saveFolder = sprintf('run%d',runNr);
    runComputationalObserverModel('conedensity', 'saveFolder', saveFolder, 'seed', runNr, 'currentFlag', true)
end

%% 2. Cone absorptions--> RGC responses
% RATIO:
% The variable ratio are integers from 1-5 that map to cone:mRGC ratios,
% assuming one mRGC = ON + OFF cell.
% Choose from: 1,2,3,4,5 --> cone:mRGC = 1:2, 1:4, 1:9, 1:16, and 1:25.
% For the manuscript, we use mRGC:cone ratio so 2:1, 1:0.5, 1:0.22,
% 1:0.125, 1:0.08. (Or 2./((1:5).^2))

% ECCEN: Eccentricity index are integers that map to different cone
% densities. Choose from 1-13 --> where 1 = highest density (at the fovea
% ~20,000 cones/deg2) and 13 is the lowest density (~40 deg eccentricity at
% nasal retina of the left eye, ~500 cones/deg2).

% 2.1 Filter cone absorptions from cone density experiment for single ratio
% and single eccentricity.
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'conedensity';
seed       = 1; % we match run1 to rng seed 1, run2 to 2.. etc.
ratio      = 2; % equal to 1:4 cone:mRGC or 1:0.5 mRGC:cone
eccen      = 5; % equal to 4.5 deg eccen
linearRGCModel(baseFolder, subFolder, expName, seed, ratio, eccen)

% 2.2 Filter cone absorptions from simple scenario with L-cone mosaic for
% single ratio and single eccentricity
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'onlyL';
expName    = 'defaultnophaseshift';
seed       = 1; % only one run, so use seed =1
ratio      = 1; % equal to 1:2 = cone:mRGC ratio, or 2:1  mRGC:cone
eccen      = 1; % equal to 4.5 deg eccen

% 2.3 Filter cone absorptions from ideal observer
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'onlyL';
expName    = 'idealobserver';
seed       = 1; % only one run, so use seed =1
ratio      = 1; % equal to 1:2 = cone:mRGC ratio, or 2:1  mRGC:cone
eccen      = 1; % equal to 4.5 deg eccen

%% 3. FIGURE 5: RGC responses --> IDEAL/SVM-Fourier/SVM-Energy performance
% For Figure 5, we demonstrate 3 classifiers:
% ideal observer, SVM-Fourier and SVM-Energy template. You can also request
% SVM-linear template percent accuracy, to compare against SVM-Energy template.

% 3.1 Ideal
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'onlyL';
expName    = 'idealobserver';
for ratio = 1:5
    getIdealObserverAccuracy(baseFolder, expName, subFolder, ratio)
end

% 3.2 SVM-Fourier
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'onlyL';
expName    = 'defaultnophaseshift';
for ratio = 1:5
    linearRGCModel_Classify(baseFolder, subFolder, expName, 1, ratio, 1);
end

% 3.3 SVM-Fourier
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'onlyL';
expName    = 'defaultnophaseshift';
for ratio = 1:5
    linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, 1, ratio, 1);
end

%% 4. FIGURE 6-7: RGC responses --> SVM-Fourier classifier performance
% For Figure 6 and 7 we show SVM-Fourier classifier performance, we
% classify separately for each eccentricity (i.e. cone density), ratio, and
% run iteration

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensity';
for eccen = 1:13
    for ratio = 1:5
        for runNr = 1:5
            subFolder = sprintf('run%d',runNr);
            seed = runNr;
            linearRGCModel_Classify(baseFolder, subFolder, expName, seed, ratio, eccen);
        end
    end
end

%% 5. Supplementary Figure  6-7: RGC responses --> SVM-Energy classifier performance
% For Suppl Figure 6 and 7 we show SVM-Energy template classifier performance, we
% classify separately for each eccentricity (i.e. cone density), ratio, and
% run iteration
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensity';
for eccen = 1:13
    for ratio = 1:5
        for runNr = 1:5
            subFolder = sprintf('run%d',runNr);
            seed = runNr;
            linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen);
        end
    end
end

%% 6. Supplementary Figure 7: Cone absorptions --> SVM-Energy template classifier performance
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
subFolder  = 'run1';
expName    = 'conedensity';
seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.

[P_svmEnergy, P_svmLinear] = absorptionsClassifierAccuracyStimTemplateWrapper...
    (baseFolder, subFolder, expName, seed)


%% 7. Combining percent accuracy across 5 runs

% 7.1 For cone absorptions only
s_combineAndMeanConeDensityExperiment

% 7.2 For cone current only
s_combineAndMeanSingleContrastCurrentAccuracyData

% 7.3 For RGC simulations classified with SVM-Fourier
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensity';
for ratio  = 1:5
    averageConeDensitySimulations(baseFolder, expName, false, ratio)
end

% 7.4 For RGC simulations classified with SVM-Energy template
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensity';
for ratio  = 1:5
    averageConeDensitySimulations(baseFolder, expName, true, ratio)
end