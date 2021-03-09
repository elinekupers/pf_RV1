function P_svm = linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen)
%
% Function to classify RGC responses computed by linearRGCmodel.m.
%
% Our classifier determines whether the responses belong to a trial with a
% 2-AFC orientation discrimination task (was Gabor clockwise or counter-
% clockwise oriented). This function applies a 2D FFT on the RGC responses
% before classifying with a linear SVM.
%
% See linearRGCModel.m to get mRGC responses.
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
% subFolder  = 'run1';
% expName    = 'conedensity';
% seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
% ratio      = 2; % cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell)
% eccen      = 5; % equal to 4.5 deg eccen

% linearRGCModel_ClassifyStimTemplate(baseFolder, subFolder, expName, seed, ratio, eccen)

%% 0. Define params
% set random number generator seed and general flags
rng(seed);

saveData = true;

% Get experimental params
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
inputType = 'absorptions';%'absorptions'; % could also be 'current'
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 1:109; % photocurrent responses are temporally delayed
end

eccentricities = expParams.eccentricities; % deg (should be 4.5 deg)

if (ratio == 5) && (any(eccen==[10,11,12,13]))
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
rgcParams.seed       = seed;

cone2RGCRatio = ratio;

% Subsampling ratio
fprintf('Ratio %d:1\n', cone2RGCRatio)

fprintf('Eccentricity %2.2f\n', eccentricities(eccen))


%% Get 2-AFC SVM Linear Classifier accuracy

% Preallocate space
P_svm = NaN(1,length(contrasts));

for c = 1:length(contrasts)
    
    % Load RGC responses
    if strcmp(expName, 'conedensity')
        load(fullfile(baseFolder, 'data',  expName, 'rgc', 'meanPoissonPadded', subFolder, sprintf('ratio%d',ratio), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c), eccentricities(eccen), inputType)), 'rgcResponse');
    else    
        inputType = 'absorptionrate';
        load(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ratio), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_%s.mat', cone2RGCRatio,  contrasts(c), inputType)), 'rgcResponse');
    end
    
    % Truncate time if needed
    if size(rgcResponse,4) > length(selectTimePoints)
        dataIn = rgcResponse(:,:,:,selectTimePoints,:);
    else
        dataIn = rgcResponse;
    end
        
    %% Get template from Gabor projected to rgc array

    stimTemplate =  getStimTemplateForSVMClassification(baseFolder, subFolder, expName, cone2RGCRatio, contrasts(c), eccentricities(eccen), selectTimePoints);
    
    % Classify!
    P_svm(c) = getClassifierAccuracyStimTemplate(dataIn, stimTemplate);
    
    fprintf('%3.2f\n',P_svm(c))
end

% Save classification results
if saveData
    if strcmp(expName,'conedensity')
        extraSubfolder = 'meanPoissonPadded';
    else
        extraSubfolder = '';
    end
    if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc',  extraSubfolder, 'stimTemplate',subFolder), 'dir');
        mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc', extraSubfolder,'stimTemplate',subFolder)); end
    parsave(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', extraSubfolder,'stimTemplate',subFolder, sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_%s_%s.mat', cone2RGCRatio, inputType, eccen, expName, subFolder)), 'P_svm',P_svm, 'rgcParams',rgcParams, 'expParams', expParams);
end

return

%% Old stuff to make template, not sure if irrelevant yet.

%     % Get template from ideal observer simulation
%     templateDir = fullfile(baseFolder, 'data', 'conedensitynonoise', 'rgc', sprintf('ratio%d',cone2RGCRatio));
%     templateFile = sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_%s.mat', cone2RGCRatio, contrasts(c), inputType);
%
%     templateNoiseless = load(fullfile(templateDir, templateFile));
%     template = templateNoiseless.rgcResponse;

%     %% Get the trials and samples (should be the data for all data sets though
%     nStimuli = size(template,5);
%     nTrials  = size(template,1) * nStimuli/2;
%     tSamples = size(template,4);
%     nrows    = size(template,2);
%     ncols    = size(template,3);
%
%     %   permute to trials x stimuli x rows x cols x time points
%     template = permute(template, [1 5 2:4]);
%
%     %   reshape to (trials x stimuli) x rows x cols x time points
%     template = reshape(template, [], nrows, ncols, tSamples);
%
%     % Label clockwise and counterclockwise trials
%     label = [ones(nTrials, 1); -ones(nTrials, 1)]; % First set is CW, second set is CCW
%
%     templateCW = template(label==1,:,:,selectTimePoints);
%     templateCW = sum(templateCW(1,:,:,:),4); % sum across all time points to have a fair comparison to the SVM results.
%     templateCCW = template(label==-1,:,:,selectTimePoints);
%     templateCCW = sum(templateCCW(1,:,:,:),4); % sum across all time points to have a fair comparison to the SVM results.
%
%     stimTemplate = squeeze(templateCW)-squeeze(templateCCW);