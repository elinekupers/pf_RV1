%% s_combineAndMeanCurrentDataConeDensityExperiment

%% 0. Set general experiment parameters
saveData   = false;
expName    = 'defaultnophaseshift'; %'conedensity';
expParams  = loadExpParams(expName, false);
[xUnits, colors, labels, coneDensity] = loadWeibullPlottingParams('current');

if strcmp(expName, 'conedensity')
    polarAngles = [0, pi/2, pi, 3*pi/2];
    polarAngleLabels = {'nasal','superior','temporal','inferior'};
else
    polarAngles = 0;
    polarAngleLabels = {'nasal'};
end

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification','current',expName);
figurePth   = fullfile(baseFolder,'figures', 'psychometricCurves', expName, 'current');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get cone ratio
lmsRatio    = expParams.cparams.spatialDensity;


P_all = [];
% First combine fovea
for pa = 1:length(polarAngles)
    
    if strcmp(expName, 'conedensity')
        fName   = sprintf('current_Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
        max(expParams.contrastLevelsPC),polarAngles(1),sprintf('%i',expParams.eyemovement'),expParams.eccentricities(5),expParams.defocusLevels,expParams.spatFreq);
    
        d = dir(fullfile(dataPth, sprintf('run*_4.5deg_%s_current',polarAngleLabels{pa})));
    
    elseif strcmp(expName, 'defaultnophaseshift')

        fName = sprintf('current_Classify_coneOutputs_contrast0.1000_pa%d_eye00_eccen%1.2f_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat', pa, expParams.eccentricities(1));
        d = dir(fullfile(dataPth, 'run*_currentAllTimePoints'));
    end
    
    
    P =[];
    for ii = 1:size(d,1)
        
        if strcmp(expName, 'defaultnophaseshift')
            tmp = load(fullfile(d(ii).folder, d(ii).name, fName));
            accuracy.P = tmp.accuracy;
        else
            accuracy = load(fullfile(d(ii).folder, d(ii).name, fName));
        end
            if size(accuracy.P,1)<size(accuracy.P,2)
            accuracy.P = accuracy.P';
        end
        
        P = [P accuracy.P];
        
    end
    
%     P_all = cat(3,P_all,P);
    P_SE = std(P,[],2)./sqrt(size(P,2));
    P_AVG = mean(P,2);
    
    figure; hold all; 
    for ii = 1:size(P,2); plot(expParams.contrastLevelsPC,P(:,ii)); end
    plot(expParams.contrastLevelsPC, P_AVG, 'k', 'lineWidth',2)
    plot(expParams.contrastLevelsPC, 80*ones(size(expParams.contrastLevelsPC)), 'k')

    if saveData
    fName   = sprintf('current_Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_AVERAGE.mat', ...
        max(expParams.contrastLevelsPC),rad2deg(polarAngles(pa)),sprintf('%i',expParams.eyemovement'),expParams.eccentricities(5),expParams.defocusLevels,expParams.spatFreq);
    
    fNameSE   = sprintf('current_Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_SE.mat', ...
        max(expParams.contrastLevelsPC),rad2deg(polarAngles(pa)),sprintf('%i',expParams.eyemovement'),expParams.eccentricities(5),expParams.defocusLevels,expParams.spatFreq);
    
        if ~exist(fullfile(dataPth, sprintf('average_4.5deg_%s_current',polarAngleLabels{pa})),'dir')
            mkdir(fullfile(dataPth, sprintf('average_4.5deg_%s_current',polarAngleLabels{pa}))); end;
        save(fullfile(dataPth,sprintf('average_4.5deg_%s_current',polarAngleLabels{pa}),fName),'P_AVG');
        save(fullfile(dataPth, sprintf('average_4.5deg_%s_current',polarAngleLabels{pa}),fNameSE),'P_SE');
    end  

    %% Bootstrap runs with replacement,
    nboot = 1000;
    bootData = bootstrp(nboot, @mean, P');
    
    % fit Weibull to each mean and get threshold
    [xUnits, ~, ~, ~, ~] = loadWeibullPlottingParams('current');
    
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
        ctrthresh(pa,ii) = ctrvar(2);
    end
    
    
    
end

if saveData
    varThresh = std(ctrthresh,[],2);
    fNameSEThresh = sprintf('varThresh_coneResponse_current_5_conedensity.mat');
    baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';

    saveFolder = fullfile(baseFolder,'data',expName,'thresholds','currentNoRGC');
    if ~exist(saveFolder, 'dir') 
        mkdir(saveFolder); 
    end
    save(fullfile(saveFolder,fNameSEThresh),'varThresh', 'ctrthresh');
end
% 
% figure(101);  clf; hold all;
% for ii = 1:size(P_all,2) 
%     subplot(411); hold all;
%     plot(expParams.contrastLevelsPC,P_all(:,ii,1));
%     set(gca,'XScale', 'log')
%     title('Nasal')
%     subplot(412); hold all;
%     plot(expParams.contrastLevelsPC,P_all(:,ii,2));
%     set(gca,'XScale', 'log')
%     title('Superior')
%     
%     subplot(413); hold all;
%     set(gca,'XScale', 'log')
%     plot(expParams.contrastLevelsPC,P_all(:,ii,3));
%     title('Temporal')
%     
%     subplot(414); hold all;
%     plot(expParams.contrastLevelsPC,P_all(:,ii,4));
%     set(gca,'XScale', 'log')
%     title('Inferior')
% end
% 
% subplot(411);
% plot(expParams.contrastLevelsPC,mean(P_all(:,:,1),2), 'k', 'lineWidth',2);
% subplot(412);
% plot(expParams.contrastLevelsPC,mean(P_all(:,:,2),2), 'k', 'lineWidth',2);
% subplot(413);
% plot(expParams.contrastLevelsPC,mean(P_all(:,:,3),2), 'k', 'lineWidth',2);
% subplot(414);
% plot(expParams.contrastLevelsPC,mean(P_all(:,:,3),2), 'k', 'lineWidth',2);
    


