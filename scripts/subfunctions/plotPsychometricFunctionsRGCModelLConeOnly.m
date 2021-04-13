function [] = plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, expName, subFolder, decisionmaker, varargin)
% Function to compute psychometric functions based on the computational
% observer model.
%
% INPUTS:
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% [subFolderName] : string defining the sub folder you want to plot from.
% [saveFig]       : boolean defining to save figures or not
% [plotAvg]      : boolean defining to plot average across experiments runs or not
%
% OUTPUTS:
% none
%
% Example 1 - ideal observer with L-cone only mosaic:
% baseFolder = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'idealobserver', 'onlyL', 'Ideal')
%
% Example 2 - SNR observer with L-cone only mosaic:
% baseFolder = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', 'onlyL', 'SNR')
%
% Example 3 - SVM observer with L-cone only mosaic:
% baseFolder = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model/';
% plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'defaultnophaseshift', 'onlyL', 'SVM')
%
% Written by EK @NYU (2020)

%% 0. Set general experiment parameters
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addRequired('subFolder', @ischar);
p.addRequired('decisionmaker', @ischar);
p.addParameter('inputType', 'absorptions', @ischar);
p.addParameter('saveFig', true, @islogical);
p.addParameter('plotAvg', false, @islogical);
p.parse(baseFolder, expName, subFolder, decisionmaker, varargin{:});

% Rename variables
baseFolder    = p.Results.baseFolder;
expName       = p.Results.expName;
subFolder     = p.Results.subFolder;
decisionmaker = p.Results.decisionmaker;
inputType     = p.Results.inputType;
saveFig       = p.Results.saveFig;
plotAvg       = p.Results.plotAvg;

% Load specific experiment parameters
expParams    = loadExpParams(expName);
[xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('rgcratios');

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification', 'rgc', 'meanPoissonPadded', decisionmaker, subFolder);
figurePth   = fullfile(baseFolder,'figures','psychometricCurves', expName, 'meanPoissonPadded', subFolder);
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
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);
fit.thresh  = 0.75;

% Set inital slope, threshold for first stage fitting
if strcmp('Ideal', decisionmaker)
    fitType = 'poly2';
    % Define a zero point (just a very small number), to plot the 0 contrast,
    % since a log-linear plot does not define 0.
    logzero = 3e-5;
    fit.init    = [5, 0.0006]; % slope, threshold at ~80%    
    scaleFactor = 1; % to put data in percentage (not needed here)
    nTotal      = 100;
elseif strcmp('SVM-Fourier', decisionmaker)     
    fitType = 'linear';
    logzero = 3e-5;
    fit.init   = [5, 0.0006]; % slope, threshold at ~80%
    scaleFactor = 1; % to put data in percentage
elseif strcmp('SVM-Energy', decisionmaker)     
    fitType = 'linear';
    logzero = 3e-5;
    fit.init   = [5, 0.0006]; % slope, threshold at ~80%
    scaleFactor = 1; % to put data in percentage
elseif strcmp('SNR', decisionmaker) 
    fitType = 'linear';
    logzero = 3e-5;
    fit.init = [3, 0.01];
    scaleFactor = 100; % to put data in percentage
end

for ratio = 1:5
    % 2. Get correct filename
    fName = sprintf('classify%s_rgcResponse_Cones2RGC%d_%s_1_%s_%s.mat', ...
        decisionmaker, ratio, inputType, expName, subFolder);

    % 3. Load performance results
    tmp = load(fullfile(dataPth,fName));
    fn = fieldnames(tmp);
    accuracy.P(ratio,:) = squeeze(tmp.(fn{1}))*scaleFactor;
end

% Transpose matrix if necessary
if size(accuracy.P,1)<size(accuracy.P,2)
    accuracy.P = accuracy.P';
end

% Resample x-axis of Weibull function to get enough low contrast values
xUnits = [linspace(min(expParams.contrastLevels),0.04, 800),linspace(0.0401, max(expParams.contrastLevels), 400)];

%% 4. Fit Weibull
count = 1;
for ii = 1:size(accuracy.P,2)
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P(:,ii), nTotal), fit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, xUnits);
    
    %% 5. Find contrast threshold
    fit.ctrthresh{count} = fit.ctrvar{count}(2);
    fit.data{count} = accuracy.P(:,ii);
    
    count = count +1;
    
end % ratios


%% 6. Visualize psychometric curves
figure(3); clf; set(gcf,'Color','w', 'Position',  [218   316   986   488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s %s', expName, decisionmaker)); hold all;

% Remove empty cells
idx = ~cellfun(@isempty, fit.ctrpred);
fit.ctrpred = fit.ctrpred(idx);

% Only plot first two, 50:50 and last 2 functions for conetypes mixed experiment
plotIdx = 1:length(fit.ctrpred);

% Loop over all functions to plot
for ii = plotIdx
    
    % What data and fit to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;
    
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', lineStyles{ii});
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
    % add zero contrast point
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
    if plotAvg
        errorbar([logzero, expParams.contrastLevels(2:end)], dataToPlot, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

set(gca, 'XScale','log','XLim',[logzero, 1], ...
         'YLim', [40 100], ...
         'TickDir','out',...
         'TickLength',[.015 .015], ...
         'XTick', [logzero, 0.0001, 0.001, 0.01, 0.1 1], ...
         'XTickLabel',sprintfc('%1.2f',[0 0.0001, 0.001, 0.01, 0.1 1]*100), ...
         'FontSize',17, 'LineWidth',2)

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels{plotIdx}, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s',expName,decisionmaker)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s.eps',expName,decisionmaker)))
    print(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_%s',expName,decisionmaker)), '-dpng')
end

% Save thresholds and fits
if ~exist(fullfile(baseFolder,'data',expName,'thresholds','meanPoissonPadded',subFolder), 'dir')
    mkdir(fullfile(baseFolder,'data',expName,'thresholds','meanPoissonPadded',subFolder)); 
end
    save(fullfile(baseFolder,'data',expName,'thresholds','meanPoissonPadded',subFolder, ...
    sprintf('cThresholds_%s_%s', decisionmaker, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');


%% 7. Plot density thresholds
if strcmp('default',expName) || strcmp('defaultnophaseshift',expName) || strcmp('idealobserver',expName)
    plotCone2RGCRatioVSThreshold(expName, fit, xThresh, 'decisionmaker', decisionmaker,'fitTypeName',fitType, 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log');
end

