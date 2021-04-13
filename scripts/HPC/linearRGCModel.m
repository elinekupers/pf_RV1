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
% To get classifier accuracy, see linearRGCModel_Classify.m
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
expParams = loadExpParams(expName);
inputType = 'absorptions'; % could be 'absorptions' or 'current'
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28; % cone absorptions are instant, therefore we truncate data to stimulus time points only 
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 1:109; % photocurrent responses are temporally delayed, so we use all time points
end
    
eccentricities = expParams.eccentricities; % deg

% We use extra (high) contrast for highest cone2rgc ratios at low cone 
% densities to reach 100% classifier accuracy.
if (ratio == 5) && (any(eccen==[10,11,12,13]))
    contrasts = [contrasts, 0.2:0.1:1];
end


% Add string to filename when using current data
if strcmp(inputType, 'current')
    preFix = 'current_';
else
    preFix = '';
end

%% Get RGC layer params

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
rgcParams.seed       = seed;
cone2RGCRatio        = ratio;

%% Compute RGC responses from linear layer

if saveData
    saveFolder = fullfile(baseFolder, 'data',  expName, 'rgc', 'meanPoissonPadded', subFolder);
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end
end

% Subsampling ratio and eccentricity
fprintf('Ratio %d:1\n', cone2RGCRatio)
fprintf('Eccentricity %2.2f\n', eccentricities(eccen))

% preallocate space
allRGCResponses = cell(1,length(contrasts));

% Load cone responses (current or absorption rates)
for c = 1:length(contrasts)
    fprintf('Contrast %1.4f\n', contrasts(c))
    
    % get filename, load cone responses
    dataFolder = fullfile(baseFolder, 'data', expName, inputType, subFolder);
    d = dir(fullfile(dataFolder,sprintf('%sOGconeOutputs_contrast%1.4f_*eccen%1.2f_*.mat', preFix, contrasts(c), eccentricities(eccen))));
    
    % Some folders are symbolic links, in that case we need an extra subfolder
    if isempty(d)
        try
            dataFolder = fullfile(baseFolder, 'data', expName, inputType, expName, subFolder);
            d = dir(fullfile(dataFolder,sprintf('%sOGconeOutputs_contrast%1.4f_*eccen%1.2f_*.mat', preFix, contrasts(c), eccentricities(eccen))));
            tmp   = load(fullfile(d.folder, d.name));
        catch ME
            rethrow(ME)
        end
    else
        tmp   = load(fullfile(d.folder, d.name));
    end
    
    fn    = fieldnames(tmp);
    if strcmp(inputType, 'absorptions')
        coneResponse = tmp.(fn{strcmpi(fn,'absorptions')});
    else
        coneResponse = tmp.(fn{strcmpi(fn,'current')}); 
    end
    coneResponse = coneResponse(:,:,:,selectTimePoints,:);
    sparams      = tmp.(fn{strcmpi(fn,'sparams')});
    
    % Define dimensions of cone mosaic and other params of stim
    sz = size(coneResponse);
    rgcParams.nTrials    = sz(1);      % #trials per stim orientation and phase
    rgcParams.cRows      = sz(2);      % #cones on y-axis
    rgcParams.cCols      = sz(3);      % #cones on x-axis
    rgcParams.timePoints = sz(4);      % ms
    rgcParams.selectTimePoints = selectTimePoints; % truncate to stimulus "on" response only
    rgcParams.fov        = sparams.sceneFOV;   % deg
    rgcParams.c          = contrasts(c);       % contrast of stim % Michelson
    rgcParams.stimSF     = expParams.spatFreq; % cpd
    
    % Do it!
    [rgcResponse, rgcarray, DoGfilter, filteredConeCurrent] = rgcLayerLinear(coneResponse, rgcParams, expParams);
    
    if saveData
        if ~exist(fullfile(saveFolder, sprintf('ratio%d',ratio)), 'dir')
            mkdir(fullfile(saveFolder, sprintf('ratio%d',ratio)));
        end
        if ~exist(fullfile(saveFolder, 'filteredOnly'), 'dir')
            mkdir(fullfile(saveFolder, 'filteredOnly'));
        end
        parsave(fullfile(saveFolder, sprintf('ratio%d',ratio), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c), eccentricities(eccen), inputType)), 'rgcResponse',rgcResponse, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
        parsave(fullfile(saveFolder, 'filteredOnly', sprintf('filteredConeCurrent_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c),eccentricities(eccen), inputType)), 'filteredConeCurrent',filteredConeCurrent, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
    end
    
    allRGCResponses{c} = rgcResponse;
end

if saveData
    parsave(fullfile(saveFolder, sprintf('rgcDoGFilter_Cones2RGC%d_%s.mat', cone2RGCRatio, inputType)), 'DoGfilter',DoGfilter, 'rgcParams', rgcParams, 'contrasts',contrasts, 'expParams',expParams);
    parsave(fullfile(saveFolder, sprintf('rgcArray_Cones2RGC%d_%s.mat', cone2RGCRatio, inputType)), 'rgcarray',rgcarray, 'rgcParams',rgcParams, 'contrasts',contrasts, 'expParams', expParams);
end


return




