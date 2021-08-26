function fH = makeFigure_ThreshVsDownsamplingLateNoiseRGCModel(...
    weibullFit, expName, expParams, figurePth, saveFig)
    
% Get downsample factors as x-axis
downsampleFactors = 2./(1:5).^2; % RGC:cone downsample ratios for 2D arrays
xticks = fliplr(downsampleFactors); % get x axis range and xtick labels
for ii = 1:length(xticks)
    xlabels_downsample{ii} = sprintf('%1.2f', xticks(ii));
end
x = downsampleFactors';

% Extract nr of eccentricities and make color guide
nrEccen     = length(expParams.eccentricities);
colors      = jet(nrEccen+1);
lateNoiseLevel = 1;

% Convert eccentricity to cone density using Curcio et al. 1990 data
conedensityLabels =  getConeDensityLabelsForPlotting(expParams);

% Set up figure
fH = figure(10); clf; set(gcf, 'Color', 'w', 'Position', [394,225,1127,580], ...
    'NumberTitle', 'off', 'Name', sprintf('Contrast threshold vs downsample factor: %s', expName));
subplot(1,3,[1 2]);
hold all

for ec = 1:nrEccen
    % Get contrast thresholds for downsampled data
    y = weibullFit.ctrthresh(5:end,ec);
    
    % Use a linear robust fit, in log-log
    f = robustfit(log10(x),log10(y));
    f_given_x = 10.^(f(2).*log10(x) + f(1));
    
    % Get R2
    R2_eccen(ec) = corr(y,f_given_x)^2;
    
    % Plot it!
    plot(x, f_given_x, 'color',colors(ec,:), 'LineWidth', 3); hold all;
    scatter(x, y, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors(ec,:), 'LineWidth',2);
    clear f gof
end

% Make plot pretty
box off; axis square; yrange = [0.001 0.01 0.1];
xlabel('Downsample factor','FontSize',20); ylabel('Contrast threshold (%)','FontSize',20)
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log', 'XTick',xticks,'XTickLabel',xlabels_downsample, ...
    'XLim', [0.06 2.3], 'YLim', [0.001 0.1], 'YTick', yrange, ...
    'YTickLabel',sprintfc('%1.1f',yrange*100), 'XGrid','on', 'YGrid','on', ...
    'XMinorGrid','on','YMinorGrid','on', 'GridAlpha',0.25,'LineWidth',0.5);
drawnow;
h = findobj(gca,'Type','line');
legend([h(end:-1:1)],conedensityLabels, 'Location','SouthWest', 'FontSize',10); legend boxoff
title(sprintf('Contrast threshold vs downsample factor'))

% Plot separate thresholds for absorptions, current, filtered
subplot(1,3,3); hold all;
preDownsampledDataToPlot = weibullFit.ctrthresh(1:3, :);

for m = 1:size(preDownsampledDataToPlot,2)
    scatter([3 2 1], preDownsampledDataToPlot(:,m), 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors(m,:), 'LineWidth',2);
end

% Make plot pretty
box off; yrange = [0.001 0.01 0.1];
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','linear', 'YScale', 'log','XTick',[1:3],'XTickLabel',fliplr({'Absorptions','Current','Filtered'}),...
    'XTickLabelRotation', 15,'XLim', [0.5 3.5], ...
    'YLim', [0.001 0.1], 'YTick',yrange,'YTickLabel',sprintfc('%1.1f',yrange*100), ...
    'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;

% Save figure if requested
if saveFig
    savefName = sprintf('ContrastThreshold_vs_Downsampling_%s_noiselevel%1.2f',expName, lateNoiseLevel);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end