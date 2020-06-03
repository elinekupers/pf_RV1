function [] = plotPsychometricFunctionsRGCModel(baseFolder, expName, subFolder, ratio, varargin)
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
% plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'run1',2)
%
%% 0. Set general experiment parameters
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('baseFolder', @ischar);
p.addRequired('expName', @ischar);
p.addRequired('subFolder', @ischar);
p.addRequired('ratio', validScalarPosNum);
p.addParameter('inputType', 'absorptionrate', @ischar);
p.addParameter('saveFig', true, @islogical);
p.addParameter('plotAvg', false, @islogical);
p.parse(baseFolder, expName, subFolder, ratio, varargin{:});

% Rename variables
baseFolder    = p.Results.baseFolder;
expName       = p.Results.expName;
subFolder     = p.Results.subFolder;
ratio         = p.Results.ratio;
inputType     = p.Results.inputType;
saveFig       = p.Results.saveFig;
plotAvg       = p.Results.plotAvg;

% Load specific experiment parameters
expParams    = loadExpParams(expName, false);
[xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams(expName);

% Where to find data and save figures
dataPth     = fullfile(baseFolder,'data',expName, 'classification', 'rgc', subFolder);
figurePth   = fullfile(baseFolder,'figures','psychometricCurves', expName, subFolder, sprintf('ratio%d', ratio));
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
fit.init   = [2, 0.01]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Get nr of conditions
nrEccen  = length(expParams.eccentricities);

% Check if different fitting accuracy for different cone types is requested

count = 1;
for eccen = 1:nrEccen
    
    
    %% 2. Get correct filename
    if plotAvg
        %                     fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        %                         max(expParams.contrastLevels), ...
        %                         sprintf('%i',expParams.eyemovement(:,em)'), ...
        %                         expParams.eccentricities(eccen), ...
        %                         expParams.defocusLevels(df), ...
        %                         expParams.spatFreq, ...
        %                         expParams.cparams.spatialDensity(lmsIdx,2), ...
        %                         expParams.cparams.spatialDensity(lmsIdx,3), ...
        %                         expParams.cparams.spatialDensity(lmsIdx,4));
        %
        %                     fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        %                         max(expParams.contrastLevels), ...
        %                         sprintf('%i',expParams.eyemovement(:,em)'), ...
        %                         expParams.eccentricities(eccen), ...
        %                         expParams.defocusLevels(df), ...
        %                         expParams.spatFreq, ...
        %                         expParams.cparams.spatialDensity(lmsIdx,2), ...
        %                         expParams.cparams.spatialDensity(lmsIdx,3), ...
        %                         expParams.cparams.spatialDensity(lmsIdx,4));
        %
        %                     SE{count} = load(fullfile(dataPth, fNameSE));
        
    else
        fName = sprintf('classifySVM_rgcResponse_Cones2RGC%d_%s_%d_%s_%s.mat', ...
            ratio, inputType, eccen, expName, subFolder);
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


%% 6. Visualize psychometric curves

figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s ratio %d', expName, ratio)); hold all;

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
set(gca, 'XTick', [logzero, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevels(2:2:end)]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels{plotIdx}, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d',expName,ratio)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d.eps',expName,ratio)))
    print(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_ratio%d',expName,ratio)), '-dpng')
end

%% 7. Plot density thresholds
if strcmp('conedensity',expName) || strcmp('eccbasedcoverage',expName)
    
    plotConeDensityVSThreshold(expName, fit, xThresh, 'fitTypeName','linear', 'saveFig', saveFig, 'figurePth', figurePth);    
end