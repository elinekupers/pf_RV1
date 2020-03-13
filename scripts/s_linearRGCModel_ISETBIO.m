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
expName   = 'default';
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
contrasts = expParams.contrastLevelsPC;
inputType = 'absorptionrate'; % could also be 'current'
eccentricity = 4.5; % deg

if strcmp(inputType, 'current')
    preFix = 'current_';
else
    preFix = '';
end

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
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = 10.1/3;             % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10)
rgcParams.DoG.wc     = 0.53;               % DoG center Gauss weight. Range: [0,1].
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = inputType;          % are we dealing with cone absorptions or current?

%% Compute RGC responses from linear layer

% preallocate space for all conditions
allRGCResponses = cell(length(cone2RGCRatios), length(contrasts));

for ii = 1:length(cone2RGCRatios)
    
    % preallocate space for current ratio
    thisRatioRGCResponses = cell(1, length(contrasts));

    % Subsampling ratio
    rgcParams.cone2RGCRatio = cone2RGCRatios(ii);
    
    % Load cone responses (current or absorption rates)
    for c = [1:10:length(contrasts), length(contrasts)]
        
        % get filename, load cone responses
        d = dir(fullfile(baseFolder, 'data', expName, subFolder,sprintf('%sOGconeOutputs_contrast%1.4f_*.mat', preFix, contrasts(c))));
        tmp   = load(fullfile(d.folder, d.name));
        fn    = fieldnames(tmp);
        coneResponse = tmp.(fn{strcmpi(fn,'absorptions')});
        sparams      = tmp.(fn{strcmpi(fn,'sparams')});
        
        % Define dimensions of cone mosaic and other params of stim
        sz = size(coneResponse);
        rgcParams.nTrials    = sz(1);      % #trials per stim orientation and phase
        rgcParams.cRows      = sz(2);      % #cones on y-axis
        rgcParams.cCols      = sz(3);      % #cones on x-axis
        rgcParams.timePoints = sz(4);      % ms
        
        rgcParams.fov        = sparams.sceneFOV;   % deg
        rgcParams.c          = contrasts(c);       % contrast of stim % Michelson 
        rgcParams.stimSF     = expParams.spatFreq; % cpd
        
        % Plot cone current if requested
        if expParams.verbose
            % Take one trial, one orientation, one phase, and mean across time points
            img = squeeze(mean(coneResponse(1,:,:,:,1),4));
            if strcmp(inputType, 'current')
                clims = [min(img(:)),0];
            else
                clims = [0,max(img(:))];
            end
            % Plot cone current retinal image
            figure(1); clf; imagesc(img);
            axis square; colormap gray; set(gca, 'CLim', clims, 'TickDir', 'out')
            title(sprintf('Cone %s, 1 trial, average across time points', inputType)); colorbar;
            xlabel('X cones'); ylabel('Y cones');
        end
        
        % Do it!
        [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(coneResponse, rgcParams);
        
        % Accumulate RGC responses
        thisRatioRGCResponses{c} = rgcResponse;
        allRGCResponses{ii,c} = rgcResponse;
    end
    
    % Accumulate DoG filters
    allFilters{ii} = DoGfilter;
    allArrays{ii} = rgcarray;
    
    if saveData
        if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc'), 'dir'); mkdir(fullfile(baseFolder, 'data',  expName, 'rgc')); end
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d_%s.mat', ii, inputType)), 'rgcResponse', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_%s.mat', ii, inputType)), 'DoGfilter', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcArray_Cones2RGC%d_%s.mat', ii, inputType)), 'rgcarray', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
    end
end





%% Get 2-AFC SVM Linear Classifier accuracy

% preallocate space
P = NaN(length(cone2RGCRatios),length(contrasts));

for ii = 1:length(cone2RGCRatios)
    
    for c = 1:size(allRGCResponses,2)
        
        % get RGC responses one ratio at a time
        dataIn = allRGCResponses{ii,c};
        
        % Save percent correct
        P(ii,c) = getClassifierAccuracy(dataIn);
        
    end
    
    if saveData
        save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classify_rgcResponse_Cones2RGC%d_%s.mat', ii, inputType)), 'P', 'rgcParams', 'expParams', '-v7.3')
    end
end

%% Plot accuracy

% load accuracy for cone current and absorptions
rgcP = load(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classify_rgcResponse_Cones2RGC1_%s.mat', inputType)));

currentClassifyData = dir(fullfile(baseFolder, 'data', expName, 'classification', 'conecurrentRV1', sprintf('current_*.mat')));
currentP = load(fullfile(currentClassifyData.folder, currentClassifyData.name));

absorptionClassifyData = dir(fullfile(baseFolder, 'data', expName, 'classification', 'conecurrentRV1', sprintf('Classify*.mat')));
absorptionP = load(fullfile(absorptionClassifyData.folder, absorptionClassifyData.name));

figure(2); clf; hold all;
plot(rgcP.expParams.contrastLevelsPC, rgcP.P, 'bo-');
plot(currentP.expParams.contrastLevelsPC, currentP.accuracy, 'ko-', 'lineWidth',2);
plot(currentP.expParams.contrastLevelsPC, absorptionP.accuracy, 'ro-', 'lineWidth',2);
xlabel('Stimulus contrast (%)'); ylabel('Accuracy (% correct)');
title('2AFC SVM classification performance')
set(gca, 'XScale', 'log', 'TickDir', 'out', 'FontSize', 16);

legend({'RGC','Cone current', 'Cone absorptions'}, 'Location', 'Best'); legend boxoff

if saveFigs
    hgexport(gcf, fullfile(pfRV1rootPath, 'figures', 'Performance_SVMClassifier_RGC_vs_Cones'))
end

