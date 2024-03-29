function [] = plotPsychometricFunctionsRGCModel(baseFolder, expName, subFolder, ratio, varargin)
% Function to compute psychometric functions based on the computational
% observer model.
%
% INPUTS:
% baseFolder      : path to project
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% subFolder       : string defining the sub folder you want to plot from.
% ratio           : integer from 1-5, to choose mRGC:cone ratio
% [inputType]     : choose from 'current' or 'absorptions' (default)
% [saveFig]       : boolean defining to save figures or not
% [plotAvg]       : boolean defining to plot average across experiments runs or not
% [meanPoissonPaddingFlag]   : boolean defining if you want to use results
%                       with (true, default) or without padding (false) before convolution)
% [stimTemplateFlag]   : boolean defining if you want use results from  
%                          SVM-Energy template observer (true) or SVM-Fourier (true, default).
% [fitTypeName]   : string defining what function to fit contrast threshold data 
%                   Choose from:
%                   - 'lowess-mesh': locally weighted regression (default)
%                   - 'linear': powerfunction (line in log-log)
%                   - 'linear-robust': same as linear, but with detecting 
%                   and down-weighting outliers before fitting
%                   - 'poly2': 2nd-degree polynomial (in log-log)
%
% OUTPUTS:
% none
%
% Example 1 - plot psychometric function for single run of simulating cone
% density variations, classified with SVM-Fourier:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'run1',2)
%
% Example 2 - plot psychometric function for average across 5 runs of 
% simulating cone density variations, classified with SVM-Fourier:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',1)
%
% Example 3 - plot psychometric function for average across 5 runs of 
% simulating cone density variations, classified with SVM-Energy template:
% baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',1, 'stimTemplateFlag', true)
%
% Written by EK @ NYU

%% 0. Set general experiment parameters
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addRequired('subFolder', @ischar);
p.addRequired('ratio', validScalarPosNum);
p.addParameter('inputType', 'absorptions', @ischar);
p.addParameter('saveFig', false, @islogical);
p.addParameter('plotAvg', true, @islogical);
p.addParameter('meanPoissonPaddingFlag', true, @islogical);
p.addParameter('stimTemplateFlag', false, @islogical);
p.addParameter('fitTypeName','lowess-mesh',@(x) ismember(x,{'linear','linear-robust','poly2','lowess-mesh'}));
p.parse(baseFolder, expName, subFolder, ratio, varargin{:});

% Rename variables
baseFolder    = p.Results.baseFolder;
expName       = p.Results.expName;
subFolder     = p.Results.subFolder;
ratio         = p.Results.ratio;
inputType     = p.Results.inputType;
saveFig       = p.Results.saveFig;
plotAvg       = p.Results.plotAvg;
meanPoissonPaddingFlag = p.Results.meanPoissonPaddingFlag;
stimTemplateFlag = p.Results.stimTemplateFlag;
fitTypeName   = p.Results.fitTypeName;

% Load specific experiment parameters
expParams    = loadExpParams(expName);
if strcmp(inputType, 'current')
    [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('current');
else
    [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams(expName);
end
% Get correct subfolder names
if meanPoissonPaddingFlag
    extraSubFolder = 'meanPoissonPadded';
else
    extraSubFolder = 'noPaddingBeforeConvolution';
end

if stimTemplateFlag
    preFix = 'SVM-Energy';
else
    preFix = 'SVM-Fourier';
end
extraSubFolder = [extraSubFolder '/' preFix];

% Where to find data and save thresholds figures
dataPth       = fullfile(baseFolder,'data',expName, 'classification', 'rgc', extraSubFolder, subFolder);
thresholdsDir = fullfile(baseFolder,'data',expName,'thresholds','rgc', extraSubFolder, subFolder);
figurePth     = fullfile(baseFolder,'figures','psychometricCurves', expName,'rgc', extraSubFolder, subFolder, sprintf('ratio%d', ratio));
if ~exist(figurePth, 'dir')
    mkdir(figurePth); 
end

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = expParams.nTrials;

%% 1. Set Weibull fitting parameters

% Prepare fit variables
fit = [];

% Predefine cell arrays
fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

% Set inital slope, threshold for first stage fitting
fit.init   = [2, 0.01]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Get nr of conditions
nrEccen  = length(expParams.eccentricities);

count = 1;
for eccen = 1:nrEccen
    
    
    %% 2. Get correct filename
    if plotAvg
        fName = sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_conedensity_AVERAGE.mat', ratio, inputType, eccen);
        fNameSE  = sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_conedensity_SE.mat', ratio, inputType, eccen);
        
        SE{count} = load(fullfile(dataPth, fNameSE));
        
    else
            fName = sprintf('classify%s_rgcResponse_Cones2RGC%d_%s_%d_%s_%s.mat', ...
            preFix, ratio, inputType, eccen, expName, subFolder);
    end
    
    %% 3. Load performance results
    
    % load model performance
    accuracy = load(fullfile(dataPth, fName));
    fn = fieldnames(accuracy);
    accuracy.P = squeeze(accuracy.(fn{1}));
    
    % Transpose matrix if necessary
    if size(accuracy.P,1)<size(accuracy.P,2)
        accuracy.P = accuracy.P';
    end
    % Allow for higher contrasts when using ratio 5 (1 mRGC for 25 cones)
    % for lowest simulated cone densities
    if (ratio == 5) && (any(eccen==[10,11,12,13]))
        if ~stimTemplateFlag && (length(expParams.contrastLevels)<length(accuracy.P))
            expParams.contrastLevels = [expParams.contrastLevels, 0.2:0.1:1];
        elseif stimTemplateFlag && (length(expParams.contrastLevels)<length(accuracy.P))
            accuracy.P((length(expParams.contrastLevels)+1):end) = [];
        end
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
    end
        
    %% 4. Fit Weibull
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P, nTotal), fit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, xUnits);
    
    %% 5. Find contrast threshold
    fit.ctrthresh{count} = fit.ctrvar{count}(2);
    fit.data{count} = accuracy.P;
    
    count = count +1;
    
end % eccen

%% Sometimes the psychometric function does not saturate and slope can't 
% be fitted. We check for those data per eccentricity and assume slope from 
% as the geomean of slopes from well fitted data.
for ii = 1:nrEccen
    allslopes(ii) = fit.ctrvar{ii}(1); 
    poorfits(ii) = allslopes(ii)<1;
end

assumedSlope = geomean(allslopes(~poorfits));
fit.poorFits = poorfits;

if any(poorfits)
    for pfIdx = find(poorfits)        
            
        ctrvar_tmp = fminsearch(@(x) ogFitWeibullFixedSlope(x, assumedSlope, expParams.contrastLevels, fit.data{pfIdx}, nTotal), fit.init(2));
        fit.ctrvar{pfIdx} = [assumedSlope, ctrvar_tmp];
        fit.ctrpred{pfIdx} = ogWeibull(fit.ctrvar{pfIdx}, xUnits);
    
        % Find contrast threshold
        fit.ctrthresh{pfIdx} = fit.ctrvar{pfIdx}(2);
    end
end
   
% There is one eccentricity for ratio 3 that just fails, we should not use
% this data point
if ~stimTemplateFlag && ratio == 3
    fit.ctrthresh{12} = NaN;
elseif stimTemplateFlag && ratio == 1
    fit.ctrthresh{2} = NaN;
elseif stimTemplateFlag && ratio == 3
    fit.ctrthresh{12} = NaN;
end


%% 6. Visualize psychometric curves
figure(3); clf; set(gcf,'Color','w', 'Position',  [218   316   986   488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s ratio %d', expName, ratio)); hold all;

% Remove empty cells
idx = ~cellfun(@isempty, fit.ctrpred);
fit.ctrpred = fit.ctrpred(idx);

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 3e-3;

% Only plot first two, 50:50 and last 2 functions for conetypes mixed experiment
plotIdx = 1:length(fit.ctrpred);

% Loop over all functions to plot
for ii = plotIdx
    
    if (ratio == 5) && (any(ii==[10,11,12,13]))
        if ~stimTemplateFlag && length(expParams.contrastLevels)<length(accuracy.P)
            expParams.contrastLevels = [expParams.contrastLevels, 0.2:0.1:1];
        elseif stimTemplateFlag && (length(expParams.contrastLevels)<length(SE{ii}.P_SE))
            SE{ii}.P_SE((length(expParams.contrastLevels)+1):end) = [];
        end
        xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
    elseif strcmp(expName, 'default')
         expParams    = loadExpParams(expName, false);
        [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('rgcratios');
    elseif strcmp(inputType, 'current')
        [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('current');
    else
        expParams    = loadExpParams(expName, false);
        [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams(expName);
    end
    
    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;
    
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', lineStyles{ii});
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
    
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
    if plotAvg
        errorbar([logzero, expParams.contrastLevels(2:end)], dataToPlot, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

% Make axes pretty
set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevels(2:2:end)]*100))
ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

% Add legend
h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels{plotIdx}, 'Location','bestoutside'); legend boxoff

% Save figures if requested
if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d',expName,ratio)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d.eps',expName,ratio)))
    print(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d',expName,ratio)), '-dpng')
end

% Save thresholds and fits
if ~exist(thresholdsDir, 'dir'); mkdir(thresholdsDir); end
save(fullfile(thresholdsDir, sprintf('cThresholds_ratio%d_%s', ratio, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh'); 

%% 7. Plot density thresholds
if strcmp('conedensity',expName)
    load(fullfile(thresholdsDir, sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptions_13_conedensity', ratio)), 'varThresh');
    plotConeDensityVSThreshold(expName, fit, xThresh, 'varThresh', varThresh', 'fitTypeName',fitTypeName, 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log', 'RGCflag', true);    
elseif strcmp('default',expName) || strcmp('defaultnophaseshift',expName) || strcmp('idealobserver',expName)
    plotCone2RGCRatioVSThreshold(expName, fit, xThresh, 'fitTypeName','linear', 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log');    
end

