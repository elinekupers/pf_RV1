%% s_linearRGCModel_ISETBIO
%
% Script to build an mRGC layer on top of the current
% computational observer model.
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
%                             'run1','seed',1,'currentFlag',false)

%% 0. Define params

% Base folder for data and figures
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/'; % can also be ogRootPath if rerunning model


saveData = true;
saveFigs = true;

% Get experimental params
subFolder = 'onlyL'; %'onlyL'
expName   = 'defaultnophaseshift'; %'idealobserver'; defaultnophaseshift
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
inputType = 'absorptionrate'; % could also be 'current'
if strcmp(inputType, 'absorptionrate')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 51:78; % photocurrent responses are temporally delayed
end

eccentricity = expParams.eccentricities; % deg (should be 4.5 deg)

if strcmp(inputType, 'current')
    preFix = 'current_';
else
    preFix = '';
end

% Get RGC layer params
% Take several cone:RGC density ratio's (for subsampling)
cone2RGCRatios = 1:5; % linear ratio

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = expParams.verbose; % print figures or not
rgcParams.saveFigs   = true;
rgcParams.expName    = expName;
rgcParams.subFolder  = subFolder;

% Define DoG Params
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;               % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = inputType;          % are we dealing with cone absorptions or current?

%% Compute RGC responses from linear layer

if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc'), 'dir')
    mkdir(fullfile(baseFolder, 'data',  expName, 'rgc'));
end

% preallocate space for all conditions
allRGCResponses = cell(length(cone2RGCRatios), length(contrasts));

for ii = 1:length(cone2RGCRatios)
    fprintf('Ratio %d:1\n', ii)
    
    % preallocate space for current ratio
    thisRatioRGCResponses = cell(1, length(contrasts));
    
    % Subsampling ratio
    rgcParams.cone2RGCRatio = cone2RGCRatios(ii);
    
    % Load cone responses (current or absorption rates)
    for c = 1:length(contrasts)
        fprintf('Contrast %1.4f\n', contrasts(c))
        % get filename, load cone responses
        d = dir(fullfile(baseFolder, 'data', expName, subFolder,sprintf('%sOGconeOutputs_contrast%1.4f_*.mat', preFix, contrasts(c))));
        tmp   = load(fullfile(d.folder, d.name));
        fn    = fieldnames(tmp);
        if strcmp(inputType, 'absorptionrate')
            coneResponse = tmp.(fn{strcmpi(fn,'absorptions')});
        else
            coneResponse = tmp.(fn{strcmpi(fn,'current')});
        end
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
            img = squeeze(mean(coneResponse(1,:,:,selectTimePoints,1),4));
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
        
        if saveData
            save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_%s.mat', ii,  contrasts(c), inputType)), 'rgcResponse', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
        end
        
        % Accumulate RGC responses
        thisRatioRGCResponses{c} = rgcResponse;
        allRGCResponses{ii,c} = rgcResponse;
    end
    
    % Accumulate DoG filters
    allFilters{ii} = DoGfilter;
    allArrays{ii} = rgcarray;
    
    if saveData
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_%s.mat', ii, inputType)), 'DoGfilter', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
        save(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcArray_Cones2RGC%d_%s.mat', ii, inputType)), 'rgcarray', 'rgcParams', 'contrasts', 'expParams', '-v7.3');
    end
end





%% Get 2-AFC SVM Linear Classifier accuracy


switch expName
    
    case 'idealobserver'
        % Preallocate space
        P_ideal = NaN(length(cone2RGCRatios),length(contrasts));
        
        % Loop over ratios
        for ii = 1:length(cone2RGCRatios)
            dataIn = allRGCResponses(ii,:);
            P_ideal(ii,:) = getIdealObserverAccuracy(dataIn, expName, subFolder, baseFolder, ii);
        end
        
        % Save classification results
        if saveData
            if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc'), 'dir'); mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc')); end
            save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classifyIdeal_rgcResponse_Cones2RGC%d_%s.mat', ii, inputType)), 'P_ideal', 'rgcParams', 'expParams', '-v7.3')
        end
        
    case 'defaultnophaseshift'
        % Preallocate space
        P_svm = NaN(length(cone2RGCRatios),length(contrasts));
        P_snr = NaN(length(cone2RGCRatios),length(contrasts));
        
        for ii = 1:length(cone2RGCRatios)
            for c = 1:length(contrasts)
                
                load(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_%s.mat', ii,  contrasts(c), inputType)), 'rgcResponse');
                
                % get RGC responses one ratio at a time
                dataIn = rgcResponse(:,:,:,selectTimePoints,:);
                
                % Get SVM classifier performance in percent correct
                P_svm(ii,c) = getClassifierAccuracy(dataIn);
            end
        end
        
        
        
        % Get SNR classifier performance in percent correct
        for ii = 1:length(cone2RGCRatios)
            for c = 1:length(contrasts)
                load(fullfile(baseFolder, 'data',  expName, 'rgc', sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_%s.mat', ii,  contrasts(c), inputType)), 'rgcResponse');

                dataIn = rgcResponse(:,:,:,selectTimePoints,:);
                %             dataIn = allRGCResponses(ii,:);
                P_snr(ii,c) = getSNRAccuracy(dataIn,expParams,inputType, ii, c, baseFolder);
            end
        end
        
        % Save classification results
        if saveData
            if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc'), 'dir'); mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc')); end
            save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classifySNR_rgcResponse_Cones2RGC%d_%s.mat', ii, inputType)), 'P_snr', 'rgcParams', 'expParams', '-v7.3')
            save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s.mat', ii, inputType)), 'P_svm', 'rgcParams', 'expParams', '-v7.3')
        end
end

%% Visualize outcome
lineColors = lines(size(P_svm,1));
figure(2); clf; hold all;
for ii = 1:size(P_svm,1)
    plot(4e-5, P_svm(ii,1), 'o','color', lineColors(ii,:), 'lineWidth',4);
    plot(contrasts(2:end), P_svm(ii,2:end),  'o-','color', lineColors(ii,:), 'lineWidth',4);
end

xlabel('Stimulus contrast (%)'); ylabel('Accuracy (% correct)');
title('2AFC SNR classification performance'); box off;
set(gca, 'XScale', 'log', 'TickDir', 'out', 'FontSize', 16);
set(gca, 'XLim', [4e-5 max(contrasts)], 'YLim', [40 100]);

% Plot accuracy from cone absorption classification for comparison
absorptionClassifyData = dir(fullfile(baseFolder, 'data', expName, 'classification', subFolder, sprintf('Classify*.mat')));
absorptionP = load(fullfile(absorptionClassifyData.folder, absorptionClassifyData.name));

plot(absorptionP.expParams.contrastLevels(2:end), absorptionP.accuracy(2:end), 'ko:', 'lineWidth',4);
plot(4e-5, absorptionP.accuracy(1), 'ko', 'lineWidth',4);

legend({'', 'RGC:Cones=1:1',...
    '', 'RGC:Cones=1:2', ...
    '', 'RGC:Cones=1:3', ...
    '', 'RGC:Cones=1:4', ...
    '', 'RGC:Cones=1:5', ...
    '', 'Cone absorptions'}, 'Location', 'Best'); legend boxoff;

hgexport(gcf, fullfile(pfRV1rootPath, 'figures', sprintf('SVM_ClassifierPerformance_RGC-absorptions_%s_%s',expName, subFolder)))
print(fullfile(pfRV1rootPath, 'figures', sprintf('SVM_ClassifierPerformance_RGC-absorptions_defaultnophaseshift_%s_%s',expName, subFolder)), '-dpng')


