function fH = makeFigure_WeibullLateNoiseRGCModel(dataPoints, fittedWeibull, ...
                        xUnits, expName, expParams, dataTypeLabel, saveFig, figurePth, lateNoiseLevel)

% Extract nr of eccentricities and make color guide
nrEccen     = length(expParams.eccentricities);
colors      = jet(nrEccen+1);

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-4;

% Convert eccentricity to cone density using Curcio et al. 1990 data
conedensityLabels =  getConeDensityLabelsForPlotting(expParams);

fH = figure;
set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], ...
    'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s %s', expName, dataTypeLabel));
hold all;
    
% Loop over all eccentricities to plot
for ii = 1:nrEccen

    % Get data point and fitted Weibull
    dataToPlot = squeeze(dataPoints(ii, :));
    fitToPlot  = squeeze(fittedWeibull(ii, :));

    % Plot all contrasts above zero
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', '-');
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');

    % Plot zero at an arbitrary small nr, because log x-scale ignores 0
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

% Make axes pretty
xticks = [0.001 0.01 0.1 1];
set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, xticks], 'XTickLabel',sprintfc('%1.1f',[0 xticks*100]))
ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);
title(sprintf('Classifier accuracy vs cone density: %s', dataTypeLabel))
h = findobj(gca,'Type','line');
legend([h(end:-2:2)],conedensityLabels, 'Location','bestoutside'); legend boxoff

% Save figure if requested
if saveFig
    savefName = sprintf('WeibullFit_contrastVSperformance_%s_%s_noiselevel%1.2f',expName, dataTypeLabel, lateNoiseLevel);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end




return