%% s_plotWeibullLateNoiseDownsampling

runnum = 1;
pth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/';
expName = 'defaultnophaseshiftlonly500';
subfolder = sprintf('run%d', runnum);
expParams = loadExpParams(expName);

lateNoiseLevel = 1;
loadStr = sprintf('classifierAccuracy_latenoiselevel%d_withDownsampling.mat',lateNoiseLevel);
load(fullfile(pth,'conecurrent', expName, subfolder, loadStr), 'PercentCorrect', 'expParams', 'rgcParams')

%% 1. Set Weibull fitting parameters

% Prepare fit variables
fit = [];

nrDataTypes = size(PercentCorrect,2);
colors = parula(nrDataTypes+1);

% Predefine cell arrays
fit.ctrpred = cell(nrDataTypes,1);
fit.ctrvar  = cell(nrDataTypes,1);
fit.ctrr2   = cell(nrDataTypes,1);
fit.data    = cell(nrDataTypes,1);

% Set inital slope, threshold for first stage fitting
fit.init   = [3, 0.01]; % slope, threshold at ~80%
fit.thresh = 0.75;

xUnits = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 500);
%% 2. Fit Weibull

for dt = 1:nrDataTypes

    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{dt} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, PercentCorrect(:,dt), 100), fit.init);

    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{dt} = ogWeibull(fit.ctrvar{dt}, xUnits);

    %% 3. Find contrast threshold
    fit.ctrthresh{dt} = fit.ctrvar{dt}(2);
    fit.data{dt}      = PercentCorrect(:,dt);

end

%% 3. Visualize
figure(3); clf; 
set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], ...
    'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); hold all;

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-4;

% Loop over all functions to plot
plotIdx = 1:length(fit.ctrpred); 
labels = {'absorptions', 'current', 'Filtered', 'LateNoise','DownSampled1','DownSampled2','DownSampled3' 'DownSampled4' 'DownSampled5'};

for ii = plotIdx
    
    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;

    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', '-');
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
    
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
end

xticks = [0.001 0.01 0.1 1];
set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, xticks], 'XTickLabel',sprintfc('%1.1f',[0 xticks*100]))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels{plotIdx}, 'Location','bestoutside'); legend boxoff
