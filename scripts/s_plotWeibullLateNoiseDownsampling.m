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
    xlabels_downsample{ii} = sprintf('%1.2f', xticks(ii)); 
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
set(gca, 'XTick',xticks,'XTickLabel',xlabels_downsample, 'XLim', [0.06 2.3]);

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
for jj = 2:4; xlabels_eccen{jj==[2:4]} = sprintf('10^%i',jj); end

% Get cone density as x-axis
x = cDensity;

% Plot it!
figure(3); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649]);
hold all
colors2 = parula(5+1);
selectDataTypes = 5:nrDataTypes;
lowess_span = 0.3;

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
set(gca, 'XTick',10.^[2:4],'XTickLabel',xlabels_eccen, 'XLim', 10.^[2.5 4.0]);

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

%% Get cone density and down sample factor at 4.5 deg

% Get mRGC data for different meridia. 
% Order = nasal, superior, temporal,inferior.
watson2015 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'mRGCWatsonISETBIO.mat'), ...
            'mRGCRFDensityPerDeg2', 'eccDeg');
assert([length(watson2015.eccDeg) == length(0:0.05:40)]);

% Curcio et al. 1990 (left eye, retina coords, same order as mRGC)
curcio1990 = load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'conesCurcioISETBIO.mat'), ...
            'conesCurcioIsetbio', 'eccDeg','angDeg');
[~,angIdx] = intersect(curcio1990.angDeg,[0,90,180,270]);
coneDensityDeg2PerMeridian= curcio1990.conesCurcioIsetbio(angIdx,:);

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2015.eccDeg==eccToCompute); % index

% Compute cone:RGC ratio
rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;
ratioAtIdx   = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% if meshfit expects cone2rgc ratio, take reciprocal for plotting 
% ratioAtIdx = (1./ratioAtIdx);

% Find contrast threshold data for all meridians: Nasal, Superior,temporal, inferior
predictedContrastThreshold = meshFit(log10(ratioAtIdx),log10(observedConesAtEccen));

%% Make 3D mesh plot

fH4 = figure(4); clf; set(gcf, 'Position', [782 44 881 756], 'Color', 'w'); clf;
ax = plot(meshFit,[X(:) Y(:)], Z(:));
ax(1).FaceLighting = 'gouraud';
ax(1).FaceColor = [1 1 1];
ax(2).Marker = 'none';

% Add labels
xlabel('mRGC : cone ratio', 'FontSize', 20)
ylabel('Cone density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Effect of late noise RGC model on contrast threshold')

zticks = [0.001:0.001:0.01, ...
          0.02:0.01:0.1, ...
          0.2:0.1:1];
ztick_labels = cell(size(zticks));
printTick = [0.01, 0.1, 1];
for ii = 1:length(printTick)
    ztick_labels{zticks==printTick(ii)} = sprintf('%2.0f',printTick(ii)*100);
end

% Make plot pretty
set(gca, 'ZLim',[-2.5 0],'FontSize', 20, 'LineWidth', 2, ...
    'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear', ...
    'XTick', log10(fliplr(downsampleFactors)), ...
    'XTickLabel', xlabels_downsample, 'XDir','reverse',...'YDir','reverse',...
    'YLim',[2 4],'YTick',[2:1:4],'YTickLabel',{'10^2','10^3','10^4',},...
    'ZTick', log10(zticks), 'ZTickLabel',ztick_labels,...
    'TickDir', 'out','View',[-134.4000   11.2000]);
grid on;
set(gca, 'GridAlpha', .2, 'ZMinorGrid', 'off', 'YMinorGrid', 'off', 'XMinorGrid', 'off')
set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1])
axis square; material shiny;

% Plot observed/biological variations in cone:mRGC ratios on mesh
hold all;
colorsRetina = {'r', 'b', 'g', 'k'};
zlift        = [0.01, 0.01, 0.01, 0.01]; % lift markers a tiny bit for visibility

for jj = 1:4
    scatter3(log10(ratioAtIdx(jj)),log10(observedConesAtEccen(jj)), ...
        predictedContrastThreshold(jj)+zlift(jj), 300, ...
        'MarkerFaceColor', colorsRetina{jj}, 'MarkerEdgeColor','k', ...
        'LineWidth',0.1, 'Marker', 'p')
end

% Save figure 
if saveFig
    fName = sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view1',lowess_span);
    hgexport(fH4, fullfile(figurePth, [fName '.eps']))
    savefig(fH4, fullfile(figurePth, [fName '.fig']))
    print(fH4, fullfile(figurePth, [fName '.png']), '-dpng')
end

% Ratio left / cone density right
set(gca, 'View', [-45.6000   11.2000])
set(gca, 'ydir', 'reverse')

if saveFig
    fName = sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view2',lowess_span);
    hgexport(fH4, fullfile(figurePth, [fName '.eps']))
    savefig(fH4, fullfile(figurePth, [fName '.fig']))
    print(fH4, fullfile(figurePth, [fName '.png']), '-dpng')
end

