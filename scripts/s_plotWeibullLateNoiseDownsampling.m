%% s_plotWeibullLateNoiseDownsampling

runnum = 3;
pth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/';
expName = 'defaultnophaseshiftlonly500';
subfolder = sprintf('run%d', runnum);
expParams = loadExpParams(expName);

lateNoiseLevel = 1;
loadStr = sprintf('classifierAccuracy_latenoiselevel%1.1f_withDownsampling.mat',lateNoiseLevel);
load(fullfile(pth,'conecurrent', expName, subfolder, loadStr), 'PercentCorrect', 'expParams', 'rgcParams')

%% 1. Set Weibull fitting parameters

% Prepare fit variables
weibullFit = [];

nrDataTypes = size(PercentCorrect,2);
colors = parula(nrDataTypes+1);

% Predefine cell arrays
weibullFit.ctrpred = cell(nrDataTypes,1);
weibullFit.ctrvar  = cell(nrDataTypes,1);
weibullFit.ctrr2   = cell(nrDataTypes,1);
weibullFit.data    = cell(nrDataTypes,1);

% Set inital slope, threshold for first stage fitting
weibullFit.init   = [3, 0.01]; % slope, threshold at ~80%
weibullFit.thresh = 0.75;

xUnits = logspace(log10(expParams.contrastLevels(2)),log10(max(expParams.contrastLevels)), 500);
%% 2. Fit Weibull

for dt = 1:nrDataTypes
    
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    weibullFit.ctrvar{dt} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, PercentCorrect(:,dt), 100), weibullFit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    weibullFit.ctrpred{dt} = ogWeibull(weibullFit.ctrvar{dt}, xUnits);
    
    %% 3. Find contrast threshold
    weibullFit.ctrthresh(dt) = weibullFit.ctrvar{dt}(2);
    weibullFit.data{dt}      = PercentCorrect(:,dt);
    
end

%% 3. Visualize
figure(3); clf;
set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], ...
    'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); hold all;

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-4;

% Loop over all functions to plot
plotIdx = 1:length(weibullFit.ctrpred);
labels = {'absorptions', 'current', 'Filtered', 'LateNoise','DownSampled1','DownSampled2','DownSampled3' 'DownSampled4' 'DownSampled5'};

for ii = plotIdx
    
    % What to plot?
    dataToPlot = weibullFit.data{ii};
    fitToPlot  = weibullFit.ctrpred{ii}*100;
    
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



%% Contrast thresholds vs downsampling
y = weibullFit.ctrthresh(5:end);
x = 1./(1:5).^2;
f = fit(x',y','power1'); 

figure; 
scatter(x,y,80,'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',2); 
hold on; 
plot(f)
set(gca,'XScale','log','YScale','log')


% 
% [fitResult, err] = polyfit(log10(xDownsampleFactorToFit),log10(y),1);
% [yFit, delta] = polyval(fitResult,log10(xDownsampleFactorToFit), err);
% R2 = 1 - (err.normr/norm(y - mean(y)))^2;
R2 = NaN;

xticks = fliplr(x);
for ii = 1:length(xticks) 
    xticklbls{ii} = sprintf('%1.2f', xticks(ii)); 
end % get x axis range

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649])
plot(x, f(x), 'r-', 'LineWidth', 3); hold all;
scatter(x, y, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',2);

% Make plot pretty
box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log')
xlabel('Downsample rate','FontSize',20); ylabel('Contrast threshold (%)','FontSize',25)
set(gca, 'XTick',xticks,'XTickLabel',xticklbls, 'XLim', [0.03 1.1]);

yrange = [0.01 0.1 1];
set(gca,'YLim', [0.005 1], 'YTick',yrange,'YTickLabel',sprintfc('%1.0f',yrange*100));
set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;

legend off; title(sprintf('Contrast threshold vs downsample rate - R^2: %1.2f', R2))
