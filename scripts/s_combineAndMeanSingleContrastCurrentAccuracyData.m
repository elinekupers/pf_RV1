%% s_combineAndMeanSingleContrastCurrentAccuracyData

% Set general experiment parameters
baseFolder = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model/';

expName           = 'conedensity';
stimTemplateFlag  = false;
currentFlag       = true;

expParams         = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);
polarAngles       = expParams.polarAngle;

% Define data folder
if stimTemplateFlag 
    subFolder    = 'SVM-Energy'; 
    templateName = '_svmEnergy'; % choose from '_svmEnergy' or '_svmLinear'
else 
    subFolder    = 'SVM-Fourier';
    templateName = [];
end
dataPth     = fullfile(baseFolder,'data',expName,'classification','current');

% Where to save data and figures
figurePth   = fullfile(ogRootPath,'figs', expName, 'current',['average' templateName]);
dataSavePth = fullfile(dataPth, subFolder, ['average' templateName]);


% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get nr of eccentricities and contrasts
nrEccen     = length(expParams.eccentricities);
nrContrasts = length(expParams.contrastLevelsPC);
lmsRatio    = expParams.cparams.spatialDensity;

nboot = 1000;
ctrthresh = NaN(nrEccen,nboot);

% Loop over eccentricities
for ec = 4:nrEccen
    P = NaN(nrContrasts,5);
    d = dir(fullfile(dataPth, 'run*'));
    for ii = 1:size(d,1)
        for c = 1:nrContrasts
            fName   = sprintf('current_Classify_coneOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
                expParams.contrastLevelsPC(c),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq);

            tmp = load(fullfile(d(ii).folder, d(ii).name, fName));

            if stimTemplateFlag
                fn = fieldnames(tmp);
                P(c,ii) = tmp.(fn{strcmpi(fn,['P' templateName])});
            else
                P(c,ii) = squeeze(tmp.P);
            end 
        end
    end
    
    figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevelsPC,P(:,ii)); end
    P_SE = std(P,[],2)./sqrt(size(P,2));
    
    P_AVG = mean(P,2);
    fName   = sprintf('current_Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        max(expParams.contrastLevelsPC),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNameSE   = sprintf('current_Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        max(expParams.contrastLevelsPC),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    
    if ~exist(dataSavePth,'dir'), mkdir(dataSavePth); end;
    save(fullfile(dataSavePth, fName),'P_AVG');
    save(fullfile(dataSavePth, fNameSE),'P_SE');
    
    %% Bootstrap runs with replacement,
    bootData = bootstrp(nboot, @mean, P');
    
    % fit Weibull to each mean and get threshold
    xUnits = linspace(min(expParams.contrastLevelsPC),max(expParams.contrastLevelsPC), 800);
    
    % Prepare fit variables
    fit = [];
    
    % Set inital slope, threshold for first stage fitting
    fit.init   = [2, 0.01]; % slope, threshold at ~80%
    fit.thresh = 0.75;

    for ii = 1:nboot
        %% 4. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        ctrvar = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevelsPC, bootData(ii,:)', nTotal), fit.init);
        
        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step.
        ctrpred = ogWeibull(ctrvar, xUnits);
        
        %% 5. Find contrast threshold
        ctrthresh(ec,ii) = ctrvar(2);
    end
end

varThresh = std(ctrthresh,[],2, 'omitnan');
fNameSEThresh = sprintf('varThresh_coneResponse_current_%d_conedensity.mat', ec);
if ~exist(fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC',subFolder), 'dir')
    mkdir(fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC', subFolder)); end
save(fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC',subFolder, fNameSEThresh),'varThresh', 'ctrthresh');



