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
% Example:
% baseFolder =
% subFolder = 
% plotPsychometricFunctionsRGCModelLConeOnly(baseFolder, 'idealobserver', subFolder, 'Ideal')
%
%% 0. Set general experiment parameters
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addRequired('subFolder', @ischar);
p.addRequired('decisionmaker', @ischar);
p.addParameter('inputType', 'absorptionrate', @ischar);
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
expParams    = loadExpParams(expName, false);
[xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('rgcratios');

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification', 'rgc',subFolder);
figurePth   = fullfile(baseFolder,'figures','psychometricCurves', expName, subFolder);
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

% Set inital slope, threshold for first stage fitting
if strcmp('SNR', decisionmaker)
    fit.init = [3, 0.01];
else
    fit.init   = [4, 0.005]; % slope, threshold at ~80%
end
fit.thresh = 0.75;


%% 2. Get correct filename

fName = sprintf('classify%s_rgcResponse_Cones2RGC5_%s.mat', ...
    decisionmaker, inputType);


%% 3. Load performance results

% load model performance
accuracy = load(fullfile(dataPth, fName));
fn = fieldnames(accuracy);
accuracy.P = squeeze(accuracy.(fn{1}));

% Transpose matrix if necessary
if size(accuracy.P,1)<size(accuracy.P,2)
    accuracy.P = accuracy.P';
end

if strcmp('SNR', decisionmaker) || strcmp('Ideal', decisionmaker)
    accuracy.P = accuracy.P.*100;
end

expParams.contrastLevels = accuracy.expParams.contrastLevels;
if strcmp('Ideal', decisionmaker)
    idx = find(expParams.contrastLevels == 0.004);
    expParams.contrastLevels = expParams.contrastLevels(1:idx);
    accuracy.P = accuracy.P(1:idx,:);
    fitType = 'poly2';
    % Define a zero point (just a very small number), to plot the 0 contrast,
    % since a log-linear plot does not define 0.
    logzero = 3e-5;

elseif strcmp('SVM', decisionmaker) 
    idx1 = find(expParams.contrastLevels == 0.001);
    idx2 = find(expParams.contrastLevels == 0.04);
    expParams.contrastLevels = expParams.contrastLevels(idx1:idx2);
    accuracy.P = accuracy.P(idx1:idx2,:);
    fitType = 'linear';
    logzero = 3e-4;

elseif strcmp('SNR', decisionmaker) 
    idx = find(expParams.contrastLevels == 0.001);
    expParams.contrastLevels = expParams.contrastLevels([idx:end]);
    accuracy.P = accuracy.P([idx:end],:);
    fitType = 'linear';
    logzero = 3e-4;

end

xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 800);


count = 1;

for ii = 1:size(accuracy.P,2)
    %% 4. Fit Weibull
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

set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.3f',[0 expParams.contrastLevels(2:2:end)]*100))

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
if ~exist(fullfile(baseFolder,'data',expName,'thresholds'), 'dir'); mkdir(fullfile(baseFolder,'data',expName,'thresholds')); end
save(fullfile(baseFolder,'data',expName,'thresholds', sprintf('cThresholds_%s_%s', decisionmaker, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');


%% 7. Plot density thresholds
if strcmp('default',expName) || strcmp('defaultnophaseshift',expName) || strcmp('idealobserver',expName)
    plotCone2RGCRatioVSThreshold(expName, fit, xThresh, 'decisionmaker', decisionmaker,'fitTypeName',fitType, 'saveFig', saveFig, 'figurePth', figurePth, 'yScale', 'log');
end

