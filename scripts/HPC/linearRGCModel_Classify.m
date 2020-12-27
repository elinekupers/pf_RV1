function [] = linearRGCModel_Classify(baseFolder, subFolder, expName, seed, ratio, eccen)
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
% subFolder  = 'run1'
% expName    = 'conedensity'
% seed       = 1; % run1 has rng seed 1, run2 has 2.. etc.
% ratio      = 2; % cone:mRGC ratio = 1:1 (actually 2:2, as one mRGC = ON + OFF cell) 
% eccen      = 5; % equal to 4.5 deg eccen

% linearRGCModel_Classify(baseFolder, subFolder, expName, seed, ratio, eccen)

%% 0. Define params
% set random number generator seed and general flags
rng(seed);

saveData = true;

% Get experimental params
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)
inputType = 'absorptions'; % can be 'absorptions' or 'current'
if strcmp(inputType, 'absorptions')
    contrasts = expParams.contrastLevels;
    selectTimePoints = 1:28;
elseif strcmp(inputType, 'current')
    contrasts = expParams.contrastLevelsPC; % PC stands for photocurrent
    selectTimePoints = 51:78; % photocurrent responses are temporally delayed
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

%% Compute RGC responses from linear layer

if ~exist(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder), 'dir')
    mkdir(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder));
end

% Subsampling ratio
fprintf('Ratio %d:1\n', cone2RGCRatio)

fprintf('Eccentricity %2.2f\n', eccentricities(eccen))


%% Get 2-AFC SVM Linear Classifier accuracy

switch expName
    
    case {'defaultnophaseshift', 'default','conedensity'}
        % Preallocate space
        P_svm = NaN(1,length(contrasts));
        
        for c = 1:length(contrasts)
            
            load(fullfile(baseFolder, 'data',  expName, 'rgc', subFolder, sprintf('ratio%d',ratio), sprintf('rgcResponse_Cones2RGC%d_contrast%1.4f_eccen%2.2f_%s.mat', cone2RGCRatio,  contrasts(c), eccentricities(eccen), inputType)), 'rgcResponse');
            if size(rgcResponse,4) > length(selectTimePoints) 
                dataIn = rgcResponse(:,:,:,selectTimePoints,:);
            else
                dataIn = rgcResponse;
            end
            
            % Get SVM classifier performance in percent correct
            P_svm(c) = getClassifierAccuracy(dataIn);
            
            fprintf('%3.2f\n',P_svm(c))
        end
        
        % Save classification results
        if saveData
            if ~exist(fullfile(baseFolder, 'data',  expName, 'classification','rgc', subFolder), 'dir'); 
                mkdir(fullfile(baseFolder, 'data',  expName, 'classification','rgc', subFolder)); end
            parsave(fullfile(baseFolder, 'data', expName, 'classification', 'rgc', subFolder, ...
                sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_%s_%s.mat', cone2RGCRatio, inputType, eccen, expName, subFolder)), ...
                'P_svm',P_svm, 'rgcParams',rgcParams, 'expParams', expParams);
        end
end
return