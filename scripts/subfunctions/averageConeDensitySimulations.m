function [] = averageConeDensitySimulations(baseFolder, ratio)

% 0. Set general experiment parameters
expName                  = 'conedensity';
expParams                = loadExpParams(expName, false);

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName,'classification','rgc');
averageDataPth = fullfile(dataPth,'average');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get nr of conditions
nrEccen     = length(expParams.eccentricities);

for ec = 1:nrEccen
    
    fName   = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptionrate_%d_conedensity_run*.mat', ratio, ec);
    
    P  =[];
    for ii = 1:5
        d = dir(fullfile(dataPth,sprintf('run%d',ii), fName));        
        
        load(fullfile(d.folder, d.name));
        if size(P_svm,1)<size(P_svm,2)
            P_svm =P_svm';
        end  
        P = [P P_svm]; 
    end
    
    if (ratio == 5) && (any(ec==[10,11,12,13]))
        expParams.contrastLevels = [expParams.contrastLevels, 0.2:0.1:1];
    end
    figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevels,P(:,ii)); end
    P_SE = std(P,[],2)./sqrt(size(P,2));
    
    P_AVG = mean(P,2);
    fNameAVG = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptionrate_%d_conedensity_AVERAGE.mat', ratio, ec);
    fNameSE  = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptionrate_%d_conedensity_SE.mat', ratio, ec);

    if ~exist(averageDataPth,'dir'), mkdir(averageDataPth); end;
    save(fullfile(averageDataPth,fNameAVG),'P_AVG');
    save(fullfile(averageDataPth,fNameSE),'P_SE');
    
    
    %% Bootstrap runs with replacement, 
    nboot = 1000;
    bootData = bootstrp(nboot, @mean, P');
    
    % fit Weibull to each mean and get threshold 
    [xUnits, ~, ~, ~, ~] = loadWeibullPlottingParams(expName);

    % Prepare fit variables
    fit = [];

    % Set inital slope, threshold for first stage fitting
    fit.init   = [2, 0.01]; % slope, threshold at ~80%
    fit.thresh = 0.75;
    
    if (ratio == 5) && (any(ec==[10,11,12,13]))
        if length(expParams.contrastLevels)<size(P,1)
            expParams.contrastLevels = [expParams.contrastLevels, 0.2:0.1:1];
        end
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
    end
    
    for ii = 1:nboot
        %% 4. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        ctrvar = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, bootData(ii,:)', nTotal), fit.init);

        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step.
        ctrpred = ogWeibull(ctrvar, xUnits);

        %% 5. Find contrast threshold
        ctrthresh(ec,ii) = ctrvar(2);
    end
end

varThresh = std(ctrthresh,[],2);
fNameSEThresh = sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptionrate_%d_conedensity.mat', ratio, ec);
if ~exist(fullfile(baseFolder,'data',expName,'thresholds'), 'dir'); mkdir(fullfile(baseFolder,'data',expName,'thresholds')); end
save(fullfile(baseFolder,'data',expName,'thresholds',fNameSEThresh),'varThresh', 'ctrthresh');
    

end


