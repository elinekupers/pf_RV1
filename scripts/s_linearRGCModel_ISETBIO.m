%% s_linearRGCModel_ISETBIO
%
% First attempt script to build an RGC layer on top of the current
% computational observer model. The RGC layer is based on the retina-V1
% model published by Bradley, Abrams and Geisler (2014) in JoV.
%
% Our RGC layer of the model has two stages, after getting the cone
% currents, first one is linear and the second is non-linear:
%   0.  Obtain cone current, either by running observer model, or loading
%   in saved cone currents.
%   1.  Cone current -> filtered by Difference of Gaussians (DoG)
%   2.  Filtered response --> sub sampling
%
%
% To recompute cone currents without eyemovements or Gabor stim phase
% shifts, use the following command:
% runComputationalObserverModel('defaultnophaseshift', 'saveFolder', ...
%                             'conecurrentRV1','seed',1,'currentFlag',true)

%% 0. Define params

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/'; % can also be ogRootPath if rerunning model

% Get experimental params
subFolder = 'conecurrentRV1';
expName   = 'defaultnophaseshift';
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
contrasts = expParams.contrastLevelsPC;

saveData = true;
saveFigs = true;

% Get RGC layer params
% Take several cone:RGC density ratio's (for subsampling)
cone2RGCRatios = 1:5; % linear ratio

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = expParams.verbose; % print figures or not
rgcParams.saveFigs   = true;

% Define DoG Params
rgcParams.DoG.kc     = 1;                  % Gauss center sigma. (as in Bradley et al. 2014 paper)
rgcParams.DoG.ks     = 10.1;               % Gauss surround sigma. Range: ks > kc. (as in Bradley et al. 2014 paper, but seems very large to us)
rgcParams.DoG.wc     = 0.53;               % DoG center Gauss weight. Range: [0,1].
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].


%% Compute RGC responses from linear layer
for ii = 1:length(cone2RGCRatios)

    % Subsampling ratio
    rgcParams.cone2RGCRatio = cone2RGCRatios(ii);
    
    % preallocate space
    allRGCResponses = cell(length(cone2RGCRatios), length(contrasts));
    
    % Load cone current
    for c = 1:length(contrasts)
    
        fName     = sprintf('current_OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat', contrasts(c));
        load(fullfile(baseFolder, 'data', expName, subFolder, fName));
        
        % Define dimensions of cone mosaic
        sz = size(current);
        rgcParams.nTrials    = sz(1);      % #trials per stim orientation and phase
        rgcParams.cRows      = sz(2);      % #cones on y-axis
        rgcParams.cCols      = sz(3);      % #cones on x-axis
        rgcParams.timePoints = sz(4);      % ms
        
        % Plot cone current if requested
        if expParams.verbose
            % Take one trial, one orientation, one phase, and mean across time points
            img = squeeze(mean(current(1,:,:,:,1),4));
            
            % Plot cone current retinal image
            figure(1); clf; imagesc(img);
            axis square; colormap gray; set(gca, 'CLim', [min(img(:)),0], 'TickDir', 'out')
            title('Cone current, 1 trial, average across time points'); colorbar;
            xlabel('X cones'); ylabel('Y cones');
        end
        
        
        % Do it!
        [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(current, rgcParams);
        
        % Accumulate RGC responses
        allRGCResponses{ii,c} = rgcResponse;
    end
    
    if expParams.verbose
        plotFFTDoG(DoGfilter,  rgcParams, expParams, saveFigs)
    end
    
    if saveData
        rgcResponses = allRGCResponses{ii,:};
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d.mat', ii)), 'rgcResponses', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
    end
    
end

%% Get 2-AFC SVM Linear Classifier accuracy

% preallocate space
P = NaN(length(cone2RGCRatios),length(contrasts));

for ii = 1:length(cone2RGCRatios)
    
    for c = 1:length(contrasts)
        
        % get RGC responses one ratio at a time
        dataIn = allRGCResponses{ii,c};
        
        % Save percent correct
        P(ii,c) = getClassifierAccuracy(dataIn);
        
    end
    
    if saveData
        save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classify_rgcResponse_Cones2RGC%d.mat', ii)), 'P', 'rgcParams', 'expParams', '-v7.3')
    end
end

%% Plot accuracy

% load accuracy for cone current and absorptions
rgcP = load(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classify_rgcResponse_Cones2RGC1.mat')));

currentClassifyData = dir(fullfile(baseFolder, 'data', expName, 'classification', 'conecurrentRV1', sprintf('current_*.mat')));
currentP = load(fullfile(currentClassifyData.folder, currentClassifyData.name));

absorptionP = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/defaultnophaseshift/run1/Classify_coneOutputs_contrast0.1000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat');

figure(2); clf; hold all;
% plot(rgcP.expParams.contrastLevelsPC, rgcP.P, 'bo-');
plot(currentP.expParams.contrastLevelsPC, currentP.accuracy, 'ko-', 'lineWidth',2);
plot(currentP.expParams.contrastLevels, absorptionP.accuracy, 'ro-', 'lineWidth',2);
xlabel('Stimulus contrast (%)'); ylabel('Accuracy (% correct)');
title('2AFC SVM classification performance')
set(gca, 'XScale', 'log', 'TickDir', 'out', 'FontSize', 16);

legend({'RGC','Cone current', 'Cone absorptions'}, 'Location', 'Best'); legend boxoff

if saveFigs
    hgexport(gcf, fullfile(pfRV1rootPath, 'figures', 'Performance_SVMClassifier_RGC_vs_Cones'))
end

