function [fH, RGC, current, absorptions] = ...
    makeFigure_ThreshVsConeDensityLateNoiseRGCModel(...
    weibullFit, whichfit, dataTypeLabels, expName, expParams, figurePth, saveFig)

% Get nr of data types
nrDataTypes = length(dataTypeLabels);

% Get downsample factors as x-axis
downsampleFactors = 2./(1:5).^2; % RGC:cone downsample ratios for 2D arrays

% Get downsample factors as labels
for ii = 1:length(downsampleFactors)
    downsamplelbls{ii} = sprintf('mRGC : cone = %1.2f : 1', downsampleFactors(ii)); %#ok<AGROW>
end

% Get cone density xticks + labels
for jj = 2:4; xlabels_eccen{jj==[2:4]} = sprintf('10^%i',jj); end %#ok<AGROW>

% Convert eccentricity to cone density using Curcio et al. 1990 data
[~, cDensity] =  getConeDensityLabelsForPlotting(expParams);

% Get colors and downsampled data types
colors2 = parula(5+1);
selectDataTypes = 5:nrDataTypes;
lateNoiseLevel = 1;

% Plot it!
fH = figure(11); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649], ...
    'NumberTitle', 'off', 'Name', sprintf('Contrast threshold vs cone density: %s', expName));
hold all

% Fit with meshgrid function
[RGC.X,RGC.Y] = meshgrid(log10(downsampleFactors),log10(cDensity));
RGC.Z = log10(weibullFit.ctrthresh(selectDataTypes,:))';

[RGC.meshFit, RGC.gof] = ogSurfaceFit(RGC.X, RGC.Y, RGC.Z, whichfit);

% Get R2
RGC.R2 = RGC.gof.rsquare;

% Store fit type
RGC.whichfit = whichfit;

% Extract single lines for plotting separate mRGC:cone ratio's
for dt = 1:length(selectDataTypes)
    RGC.data(dt,:) = weibullFit.ctrthresh(selectDataTypes(dt),:)';
    RGC.lineFit(dt,:) = 10.^RGC.meshFit(RGC.X(:,dt),RGC.Y(:,dt));
    plot(cDensity, RGC.lineFit(dt,:), 'color',colors2(dt,:), 'LineWidth', 3); hold all;
    scatter(cDensity, RGC.data(dt,:), 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors2(dt,:), 'LineWidth',2);
end

%% FILTERED: Fit filtered current by RGC DoGs
filtered.data = weibullFit.ctrthresh(3,:)';
[filtered.X,filtered.Y] = meshgrid(ones(11,1),log10(cDensity));
filtered.Z              = repmat(log10(filtered.data),1,11);
[filtered.meshFit, filtered.gof] = ogSurfaceFit(filtered.X, filtered.Y, filtered.Z, whichfit);

% Extract single lines for separate ratio's
filtered.lineFit = 10.^filtered.meshFit(filtered.X(:,1), filtered.Y(:,1));
filtered.R2      = filtered.gof.rsquare;

% Store fit type
filtered.whichfit = whichfit;

% Plot fit and markers
colorsGray = [0.7 0.7 0.7; 0.3 0.3 0.3];
plot(cDensity, filtered.lineFit, 'color', colorsGray(1,:), 'LineWidth', 2); hold all;
scatter(cDensity, filtered.data, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colorsGray(1,:), 'LineWidth',2);

%% CURRENT: Fit and plot current data
current.data          = weibullFit.ctrthresh(2,:)';
[current.X,current.Y] = meshgrid(ones(11,1),log10(cDensity)); % make 2D dummy grid
current.Z             = repmat(log10(current.data),1,11);
[current.meshFit, current.gof] = ogSurfaceFit(current.X, current.Y, current.Z, whichfit);

% Extract single line from dummy grid
current.lineFit = 10.^current.meshFit(current.X(:,1),current.Y(:,1));
current.R2      = current.gof.rsquare;

% Store fit type
current.whichfit = whichfit;

% Plot fit and markers
plot(cDensity, current.lineFit, 'color', colorsGray(2,:), 'LineWidth', 2); hold all;
scatter(cDensity, current.data, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colorsGray(2,:), 'LineWidth',2);

%% ABSORPTIONS: fit absorptions with meshfit
% or as straight line in log-log (i.e., powerlaw)
absorptions.data = weibullFit.ctrthresh(1,:)';
% [f_absorptions, gof1] = fit(cDensity,y_absorptions,'power1');
% R2_absorptions = gof1.rsquare; % Coefficient of Determination
% plot(cDensity, f_absorptions(cDensity), 'color','k', 'LineWidth', 3); hold all;
[absorptions.X, absorptions.Y] = meshgrid(ones(11,1),log10(cDensity));
absorptions.Z                  = repmat(log10(absorptions.data),1,11);
[absorptions.meshFit, absorptions.gof] = ogSurfaceFit(absorptions.X, absorptions.Y, absorptions.Z, whichfit);

% Extract single lines from dummy grid
absorptions.lineFit = 10.^absorptions.meshFit(absorptions.X(:,1),absorptions.Y(:,1));
absorptions.R2      = absorptions.gof.rsquare;

% Store fit type
absorptions.whichfit = whichfit;

plot(cDensity, absorptions.lineFit, 'color', 'k', 'LineWidth', 2); hold all;
scatter(cDensity, absorptions.data, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',2);

% Make plot pretty
box off; axis square; yrange = [0.001 0.01 0.1];
xlabel('Cone density (cones/deg^2)','FontSize',20); ylabel('Contrast threshold (%)','FontSize',20)
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log', 'XTick',10.^[2:4],'XTickLabel',xlabels_eccen,...
    'XLim', 10.^[2.5 4.0], 'YLim', [min(yrange) max(yrange)], 'YTick',yrange, ...
    'YTickLabel',sprintfc('%1.1f',yrange*100), ...
    'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
h = findobj(gca,'Type','line');
allDataTypeLabels = {downsamplelbls{:}, 'Filtered by RGC DoG','Cone current','Cone absorptions'};
legend([h(end:-1:1)],allDataTypeLabels, 'Location','bestoutside'); legend boxoff
title(sprintf('Contrast threshold vs cone density'))

% % Arrow at 4.5 deg eccentricity
% xArrow = [0.35 0.35]; % coordinates normalized to figure size, will change if figure size changes
% yArrow = [0.75 0.65];
% annotation('textarrow',xArrow,yArrow,'String',['4.5' char(176) ' eccen'], ...
%     'FontSize',13,'LineWidth',2)

% Save figure if requested
if saveFig
    savefName = sprintf('ContrastThreshold_vs_Conedensity_%s_noiselevel%1.2f_%s', ...
        expName, lateNoiseLevel, whichfit);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end
