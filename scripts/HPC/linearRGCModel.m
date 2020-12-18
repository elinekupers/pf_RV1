function []=linearRGCModel(baseFolder, subFolder, expName, seed, ratio, eccen)
%
% Function to build an RGC layer on top of the current
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
%
% INPUTS
% baseFolder    :   folder where rgc data live 
% subFolder     :   subfolder name when running multiple iterations.
% expName       :   string with experiment name, e.g. 'conedensity' (for
%                   all options see loadExpParams.m)
% seed          :   random number generator seed, in manuscript, we use the
%                   following rule: run1 has rng seed 1, run2 has 2.. etc.
% ratio         :   integer between 1-5 to indicate cone:mRGC ratio
%                   because one mRGC = ON + OFF cell, 1 = 2:1, 2 = 1:1, 3 = 0.67:1, 4 = 0.5:1, 5 = 0.4:1 
% eccen         :   index of eccentricity vector, for conedensity exp, eccentricies = [0 0.5 1 2 4.5 5 10:5:40] 
%
% Example:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% subFolder  = 'run1'
% expName    = 'conedensity'
% seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
% ratio      = 2; % cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell) 
% eccen      = 5; % equal to 4.5 deg eccen

% linearRGCModel(baseFolder, subFolder, expName, seed, ratio, eccen)

%% 0. Define params

% Base folder for data and figures
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/'; %'/scratch/ek99/pf_RV1';

rng(seed);

saveData = true;

% Get experimental params
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
inputType = 'current'; % could be 'absorptions' or 'current'
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 51:78; % photocurrent responses are temporally delayed
end
    
eccentricities = expParams.eccentricities; % deg

if (ratio == 5) && (any(eccen==[11,12,13]))
    contrasts = [contrasts, 0.2:0.1:1];
end

if strcmp(inputType, 'current')
    preFix = 'current_';
else
    preFix = '';
end

% Get RGC layer params

% Define general RGC params
rgcParams = struct();
rgcParams.verbose    = false; % print figures or not
rgcParams.saveFigs   = true;
rgcParams.expName    = expName;
rgcParams.subFolder  = subFolder;

% Define DoG Params
rgcParams.DoG.kc     = 1/3;                % Gauss center sigma. (Bradley et al. 2014 paper has kc =1)
rgcParams.DoG.ks     = rgcParams.DoG.kc*6;  % Gauss surround sigma. Range: ks > kc. (Bradley et al. 2014 paper has ks = 10.1)
rgcParams.DoG.wc     = 0.64;               % DoG center Gauss weight. Range: [0,1]. (Bradley et al. 2014 paper has ws = 0.53)
rgcParams.DoG.ws     = 1-rgcParams.DoG.wc; % DoG surround Gauss weight. Range: [0,1].
rgcParams.inputType  = inputType;          % are we dealing with cone absorptions or current?

rgcParams.cone2RGCRatio = ratio;           % linear ratio
rgcParams.seed          = seed;
cone2RGCRatio           = ratio;
%% Compute RGC responses from linear layer

if saveData
    if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder), 'dir')
        mkdir(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder));
    end
end

% Subsampling ratio
fprintf('Ratio %d:1\n', cone2RGCRatio)

fprintf('Eccentricity %2.2f\n', eccentricities(eccen))

allRGCResponses = cell(1,length(contrasts));

% Load cone responses (current or absorption rates)
for c = 1:length(contrasts)
    fprintf('Contrast %1.4f\n', contrasts(c))
    % get filename, load cone responses
    d     = dir(fullfile(baseFolder, 'data', expName, inputType, expName, subFolder,sprintf('%sOGconeOutputs_contrast%1.3f_*eccen%1.2f_*.mat', preFix, contrasts(c), eccentricities(eccen))));
    tmp   = load(fullfile(d.folder, d.name));
    fn    = fieldnames(tmp);
    if strcmp(inputType, 'absorptions')
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
    
    % Do it!
    [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(coneResponse, rgcParams);
    
    if saveData
        if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ratio)), 'dir')
            mkdir(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ratio)));
        end
        if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, 'filteredOnly'), 'dir')
            mkdir(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, 'filteredOnly'));
        end
        parsave(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ratio), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c), eccentricities(eccen), inputType)), 'rgcResponse',rgcResponse, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
        parsave(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, 'filteredOnly', sprintf('filteredConeCurrent_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c),eccentricities(eccen), inputType)), 'filteredConeCurrent',filteredConeCurrent, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
    end
    
    allRGCResponses{c} = rgcResponse;
end

if saveData
    parsave(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('rgcDoGFilter_Cones2RGC%d_%s.mat', cone2RGCRatio, inputType)), 'DoGfilter',DoGfilter, 'rgcParams', rgcParams, 'contrasts',contrasts, 'expParams',expParams);
    parsave(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('rgcArray_Cones2RGC%d_%s.mat', cone2RGCRatio, inputType)), 'rgcarray',rgcarray, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
end



% 
% 
% %% Get 2-AFC SVM Linear Classifier accuracy
% 
% switch expName
%     
%     case 'idealobserver'
%         % Preallocate space
%         P_ideal = NaN(1,length(contrasts));
%         
%         % Loop over ratios
%         dataIn = allRGCResponses{:};
%         P_ideal = getIdealObserverAccuracy(dataIn, expName, subFolder, baseFolder, allRGCResponses);
%         
%         % Save classification results
%         if saveData
%             if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc'), 'dir'); mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc')); end
%             save(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', sprintf('classifyIdeal_rgcResponse_Cones2RGC%d_%s.mat', cone2RGCRatio, inputType)), 'P_ideal', 'rgcParams', 'expParams', '-v7.3')
%         end
%         
%     case {'defaultnophaseshift', 'default','conedensity'}
%         % Preallocate space
%         P_svm = NaN(1,length(contrasts));
%         
%         for c = 1:length(contrasts)
%             dataIn = allRGCResponses{c};
%             
%             % Get SVM classifier performance in percent correct
%             P_svm(c) = getClassifierAccuracy(dataIn(:,:,:,selectTimePoints,:));
%             fprintf('%s\n',P_svm(c))
%         end
%         
%         % Save classification results
%         if saveData
%             if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc', subFolder), 'dir'); mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc', subFolder)); end
%             parsave(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', subFolder, sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_%s_%s.mat', cone2RGCRatio, inputType, eccen, expName, subFolder)), 'P_svm',P_svm, 'rgcParams',rgcParams, 'expParams', expParams);
%         end
% end

return




