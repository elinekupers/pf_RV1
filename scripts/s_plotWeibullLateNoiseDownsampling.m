%% s_plotWeibullLateNoiseDownsampling

%% 0. Set general parameters
runnum            = 3;
pth               = '/Volumes/server/Projects/PerformanceFieldsIsetBio/data/';
expName           = 'conedensitynophaseshiftlonly500';
subfolder         = sprintf('run%d', runnum);
expParams         = loadExpParams(expName);
saveFig           = true;
lateNoiseLevel    = 1;
dataTypeLabels    = {'absorptions', 'current', 'Filtered', 'LateNoise',...
                    'DownSampled1','DownSampled2','DownSampled3' 'DownSampled4' 'DownSampled5'};
[conedensityLabels, cDensity] =  getConeDensityLabelsForPlotting(expParams);

if saveFig
    figurePth  = fullfile('/Volumes/server/Projects/PerformanceFields_RetinaV1Model/', ...
            'figures','psychometricCurves', expName, 'current', subfolder);
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
end

%% 1. Set Weibull parameters
% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-4;

nrDataTypes = length(dataTypeLabels);
nrEccen     = length(expParams.eccentricities);
colors      = parula(nrEccen+1);
xUnits      = logspace(log10(logzero),log10(max(expParams.contrastLevels)), 500);

% Predefine cell arrays
weibullFit = struct();
weibullFit.ctrvar       = NaN(nrDataTypes,nrEccen, 2); % estimated variables
weibullFit.ctrpred      = NaN(nrDataTypes,nrEccen, length(xUnits)); % estimated fine sampled Weibull prediction
weibullFit.data         = NaN(nrDataTypes,nrEccen, length(expParams.contrastLevels));
weibullFit.ctrthresh    = NaN(nrDataTypes,nrEccen, 1); % final threshold

% Set inital slope, threshold for first stage fitting
weibullFit.init   = [3, 0.01]; % slope, threshold at ~80%
weibullFit.thresh = 0.75;

%% 3. Fit!
for eccen = 1:nrEccen
    loadStr = sprintf('classifierAccuracy_latenoiselevel%1.1f_eccen%1.2f_withDownsampling.mat', ...
                    lateNoiseLevel, expParams.eccentricities(eccen));
    load(fullfile(pth,'conecurrent', expName, subfolder, loadStr), 'PercentCorrect', 'expParams', 'rgcParams')
    
    %% 2. Fit Weibull
    for dt = 1:nrDataTypes

        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        weibullFit.ctrvar(dt, eccen, :) = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, ...
            PercentCorrect(:,dt), 100), weibullFit.init);

        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step, and finely sampled x-axis 
        weibullFit.ctrpred(dt, eccen,:) = 100*ogWeibull(weibullFit.ctrvar(dt, eccen,[1,2]), xUnits); % multiply by 100 to get percent correct

        %% 3. Find contrast threshold
        weibullFit.ctrthresh(dt, eccen) = weibullFit.ctrvar(dt, eccen, 2);
        weibullFit.data(dt, eccen, :)   = PercentCorrect(:,dt);

    end
end

%% 4. Visualize weibulls per cone density and datatype

for dt = 1:nrDataTypes

    figure;
    set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], ...
        'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); 
    hold all;

    % Loop over all eccentricities to plot
    for ii = 1:nrEccen

        % Get data point and fitted Weibull
        dataToPlot = squeeze(weibullFit.data(dt, ii, :));
        fitToPlot  = squeeze(weibullFit.ctrpred(dt, ii, :)); 
        
        % Plot all contrasts above zero
        plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', '-');
        scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:), 'filled');
        
        % Plot zero at an arbitrary small nr, because log x-scale ignores 0
        plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))

    end

    xticks = [0.001 0.01 0.1 1];
    set(gca, 'XScale','log','XLim',[logzero, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
    set(gca, 'XTick', [logzero, xticks], 'XTickLabel',sprintfc('%1.1f',[0 xticks*100]))

    ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
    xlabel('Stimulus Contrast (%)', 'FontSize',17);
    title(sprintf('Classifier accuracy vs cone density: %s', dataTypeLabels{dt}))
    h = findobj(gca,'Type','line');
    legend([h(end:-2:2)],conedensityLabels, 'Location','bestoutside'); legend boxoff
    
    if saveFig
        savefName = sprintf('WeibullFit_contrastVSperformance_%s_%s_noiselevel%1.2f',expName, dataTypeLabels{dt}, lateNoiseLevel);
        savefig(fullfile(figurePth,[savefName '.fig']))
        hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
        print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
    end
end
    
%% Plot contrast thresholds vs downsampling

% Get downsample factors as x-axis
downsampleFactors = 2./(1:5).^2; 
xticks = fliplr(downsampleFactors); % get x axis range and xtick labels
for ii = 1:length(xticks) 
    xticklbls{ii} = sprintf('%1.2f', xticks(ii)); 
end 
x = downsampleFactors';

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649]);
hold all

for ec = 1:nrEccen
    y = weibullFit.ctrthresh(5:end,ec);
    
    % Use a linear robust fit, in log-log
    f = robustfit(log10(x),log10(y));
    f_given_x = 10.^(f(2).*log10(x) + f(1));
    
    % Get R2
    R2_eccen(ec) = corr(y,f_given_x)^2;
    
%     [f, gof] = fit(x,y,'power1', 'Robust','Bisquare'); 
%     R2_eccen(ec) = gof.rsquare; % Coefficient of Determination

%     plot(x, f(x), 'color',colors(ec,:), 'LineWidth', 3); hold all;
    plot(x, f_given_x, 'color',colors(ec,:), 'LineWidth', 3); hold all;
    scatter(x, y, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors(ec,:), 'LineWidth',2);
    clear f gof
end

% Make plot pretty
box off; axis square;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log')
xlabel('Downsample factor','FontSize',20); ylabel('Contrast threshold (%)','FontSize',20)
set(gca, 'XTick',xticks,'XTickLabel',xticklbls, 'XLim', [0.06 2.3]);

yrange = [0.01 0.1 1];
set(gca,'YLim', [0.001 1], 'YTick',yrange,'YTickLabel',sprintfc('%1.0f',yrange*100));
set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
h = findobj(gca,'Type','line');
legend([h(end:-1:1)],conedensityLabels, 'Location','bestoutside'); legend boxoff
title(sprintf('Contrast threshold vs downsample factor'))

if saveFig
    savefName = sprintf('ContrastThreshold_vs_Downsampling_%s_noiselevel%1.2f',expName, lateNoiseLevel);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end

%% Plot contrast thresholds vs cone density

% Get downsample factors as labels
for ii = 1:length(downsampleFactors) 
    downsamplelbls{ii} = sprintf('mRGC : cone = %1.2f : 1', downsampleFactors(ii)); 
end

% Get cone density xticks + labels
for jj = 2:4; xlabels{jj==[2:4]} = sprintf('10^%i',jj); end

% Get cone density as x-axis
x = cDensity;

% Plot it!
figure(3); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649]);
hold all
colors2 = parula(5+1);
selectDataTypes = 5:nrDataTypes;
lowess_span = 0.15;

% Fit with meshgrid function (using dummy 2D grid)
[X,Y] = meshgrid(log10(downsampleFactors),log10(x));
Z = log10(weibullFit.ctrthresh(selectDataTypes,:))';
[meshFit, gof] = fit([X(:) Y(:)], Z(:), 'lowess','span',lowess_span);

% Get R2
R2_dt = gof.rsquare;

for dt = 1:length(selectDataTypes)    
    
    y = weibullFit.ctrthresh(selectDataTypes(dt),:)';
%     [f, gof] = fit(x,y,'power1'); 
%     R2_dt(ec) = gof.rsquare;
    
%     plot(x, f(x), 'color',colors2(dt,:), 'LineWidth', 3); hold all;
%     scatter(x, y, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors2(dt,:), 'LineWidth',2);

    % Extract single lines for separate ratio's
    f_given_x = 10.^meshFit(X(:,dt),Y(:,dt));
    plot(x, f_given_x, 'color',colors2(dt,:), 'LineWidth', 3); hold all;
    scatter(x, y, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors2(dt,:), 'LineWidth',2);

end

% Make plot pretty
box off; axis square;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log')
xlabel('Cone density (cones/deg^2)','FontSize',20); ylabel('Contrast threshold (%)','FontSize',20)
set(gca, 'XTick',10.^[2:4],'XTickLabel',xlabels, 'XLim', 10.^[2.5 4.0]);

yrange = [0.01 0.1 1];
set(gca,'YLim', [0.001 1], 'YTick',yrange,'YTickLabel',sprintfc('%1.0f',yrange*100));
set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
h = findobj(gca,'Type','line');
legend([h(end:-1:1)],downsamplelbls, 'Location','bestoutside'); legend boxoff
title(sprintf('Contrast threshold vs cone density'))

% Arrow at 4.5 deg eccentricity
xArrow = [0.345 0.345];
yArrow = [0.75 0.65];
annotation('textarrow',xArrow,yArrow,'String',['4.5' char(176) ' eccen'], ...
            'FontSize',13,'LineWidth',2)

if saveFig
    savefName = sprintf('ContrastThreshold_vs_Conedensity_%s_noiselevel%1.2f_lowessfit%1.1f', ...
        expName, lateNoiseLevel, lowess_span);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end