function [] = averageConeDensitySimulations(baseFolder, expName, stimTemplateFlag, ratio)

% 0. Set general experiment parameters
expParams                = loadExpParams(expName, false);

% Where to find data and save figures
subFolder   = 'meanPoissonPadded';
if stimTemplateFlag
    subFolder = [subFolder '/SVM-Energy']; 
    templateName = '_svmLinear'; % choose from '_svmEnergy' or '_svmLinear'
else 
    templateName = '/SVM-Fourier';
end

dataPth     = fullfile(baseFolder,'data',expName,'classification','rgc',subFolder);
averageDataPth = fullfile(dataPth,['average' templateName]);

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get nr of conditions
nrEccen     = length(expParams.eccentricities);

if strcmp(expName, 'conedensity')
    extraSubfolder = 'run';
elseif strcmp(expName, 'defaultnophaseshift')
    extraSubfolder = 'onlyL';
end

for ec = 1:nrEccen
    

    fName   = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_%s_%s*.mat', ratio, ec, expName, extraSubfolder);
    
    P = [];
    % Loop over runs (i.e. experiment iterations)
    for ii = 1:5
        
        if strcmp(expName, 'conedensity')
            extraSubfolder = sprintf('run%d',ii);
        elseif strcmp(expName, 'defaultnophaseshift')
            extraSubfolder = 'onlyL';
        end
   
        % Get file and load it
        d = dir(fullfile(dataPth, extraSubfolder, fName));        
        tmp = load(fullfile(d.folder, d.name));
        
        % Get percent correct  
        if stimTemplateFlag
            fn = fieldnames(tmp);
            P_svm = tmp.(fn{strcmpi(fn,['P' templateName])});
        end 
        
        % Transpose if necessary
        if size(P_svm,1)<size(P_svm,2)
            P_svm =P_svm';
        end
        
        % Accumulate
        P = [P P_svm];
    end
    
    % Deal with higher contrast levels for high cone to RGC ratio and low
    % cone densities
    if ~stimTemplateFlag && (ratio == 5) && (any(ec==[10,11,12,13]))
        if length(expParams.contrastLevels) < size(P,1)
            expParams.contrastLevels = [expParams.contrastLevels, 0.2:0.1:1];
        end
    elseif stimTemplateFlag && (ratio == 5) && (any(ec==[10,11,12,13]))
        if length(expParams.contrastLevels) < size(P,1)
            P = P(1:length(expParams.contrastLevels),:);
        else
            P = P(1:length(setdiff(expParams.contrastLevels,0.2:0.1:1)),:);
        end
    end
    
    % Plot accuracy
    figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevels,P(:,ii)); end
    
    % Compute standard error and mean
    P_SE = std(P,[],2)./sqrt(size(P,2));
    P_AVG = mean(P,2);
    
    % Save sample mean output
    fNameAVG = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_AVERAGE.mat', ratio, ec);
    fNameSE  = sprintf('classifySVM_rgcResponse_Cones2RGC%d_absorptions_%d_conedensity_SE.mat', ratio, ec);
    if ~exist(averageDataPth,'dir'), mkdir(averageDataPth); end;
    save(fullfile(averageDataPth,fNameAVG),'P_AVG');
    save(fullfile(averageDataPth,fNameSE),'P_SE');
    
    %% Bootstrap runs with replacement, 
    nboot = 1000;
    bootData{ec} = bootstrp(nboot, @mean, P');
    
    % fit Weibull to each mean and get threshold 
    [xUnits, ~, ~, ~, ~] = loadWeibullPlottingParams(expName);

    % Prepare fit variables
    fit = [];

    % Set inital slope, threshold for first stage fitting
    fit.init   = [2, 0.01]; % slope, threshold at ~80%
    fit.thresh = 0.75;
    
    if (ratio == 5) && (any(ec==[10,11,12,13]))
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
    end
    
    for boot = 1:nboot
        %% 4. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        ctrvar(ec,:) = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, bootData{ec}(boot,:)', nTotal), fit.init);

        %% 5. Find contrast threshold
        ctrthresh(ec,boot) = ctrvar(ec,2);
    end
end

% check for unfitted data per eccentricity 
poorfits = ctrvar(:,1)<1;
assumedSlope = geomean(ctrvar(~poorfits,1));

if any(poorfits)
    for pfIdx = find(poorfits')        
        for boot = 1:nboot
            ctrvar_tmp = fminsearch(@(x) ogFitWeibullFixedSlope(x, assumedSlope, expParams.contrastLevels, bootData{pfIdx}(boot,:)', nTotal), fit.init(2));
            ctrthresh(pfIdx,boot) = ctrvar_tmp;
        end
    end
end


varThresh = std(ctrthresh,[],2);
fNameSEThresh = sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptions_%d_%s.mat', ratio, ec, expName);
thresholdDir = fullfile(baseFolder,'data',expName,'thresholds',subFolder, ['average' templateName]);
if ~exist(thresholdDir, 'dir'); mkdir(thresholdDir); end
save(fullfile(thresholdDir,fNameSEThresh),'varThresh', 'ctrthresh');
    

end


