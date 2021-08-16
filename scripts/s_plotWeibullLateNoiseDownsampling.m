%% s_plotFiguresLateNoiseRGCmodel

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
downsampleFactors = 2./(1:5).^2; 
expParams.eccentricities = [1 2 3 4 4.5 5 6 7 10:5:40]; % deg
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
colors      = jet(nrEccen+1);
xUnits      = logspace(log10(logzero),log10(max(expParams.contrastLevels)), 500);

% Predefine cell arrays
weibullFit = struct();
weibullFit.ctrvar       = NaN(nrDataTypes,nrEccen, 2); % estimated variables
weibullFit.ctrpred      = NaN(nrDataTypes,nrEccen, length(xUnits)); % estimated fine sampled Weibull prediction
weibullFit.data         = NaN(nrDataTypes,nrEccen, length(expParams.contrastLevels));
weibullFit.ctrthresh    = NaN(nrDataTypes,nrEccen, 1); % final threshold
weibullFit.init = [3, 0.01];


%% 2. Fit!
for eccen = 1:nrEccen
    loadStr = sprintf('classifierAccuracy_latenoiselevel%1.1f_eccen%1.2f_withDownsampling.mat', ...
                    lateNoiseLevel, expParams.eccentricities(eccen));
    load(fullfile(pth,'conecurrent', expName, subfolder, loadStr), 'PercentCorrect')
    
    % Fit Weibull
    for dt = 1:nrDataTypes
%         if dt == 1 
%             weibullFit.init = [2, 0.001];
%         end
        accuracy = PercentCorrect(:,dt);
        w = fitWeibullToClassifierAccuracy(expParams.contrastLevels, accuracy, weibullFit.init, xUnits);
        
        weibullFit.ctrvar(dt,eccen, :) = w.ctrvar;
        weibullFit.ctrpred(dt, eccen, :) = w.ctrpred;
        weibullFit.data(dt, eccen,:) = w.data;
        weibullFit.ctrthresh(dt,eccen,:) = w.ctrthresh;
        
    end
end

%% 3. Visualize weibulls per cone density and datatype

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
    
%% 4. Plot contrast thresholds vs downsampling

% Get downsample factors as x-axis
xticks = fliplr(downsampleFactors); % get x axis range and xtick labels
for ii = 1:length(xticks) 
    xlabels_downsample{ii} = sprintf('%1.2f', xticks(ii)); 
end 
x = downsampleFactors';

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [394,225,1127,580]);
subplot(1,3,[1 2]);
hold all

for ec = 1:nrEccen
    y = weibullFit.ctrthresh(5:end,ec);
    
    % Use a linear robust fit, in log-log
    f = robustfit(log10(x),log10(y));
    f_given_x = 10.^(f(2).*log10(x) + f(1));
    
    % Get R2
    R2_eccen(ec) = corr(y,f_given_x)^2;
    
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
legend([h(end:-1:1)],conedensityLabels, 'Location','SouthWest', 'FontSize',10); legend boxoff
title(sprintf('Contrast threshold vs downsample factor'))

% Plot separate points for absorptions, current, filtered
subplot(1,3,3); hold all;
preDownsampledDataToPlot = weibullFit.ctrthresh(1:3, :);

for m = 1:size(preDownsampledDataToPlot,2)
    scatter([1:3], preDownsampledDataToPlot(:,m), 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colors(m,:), 'LineWidth',2);
end


% Make plot pretty
box off; 
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','linear', 'YScale', 'log')
xlabel('Data types','FontSize',20);
set(gca, 'XTick',[1:3],'XTickLabel',{'Absorptions','Current','Filtered'},...
    'XTickLabelRotation', 15);
set(gca,'XLim', [0.5 3.5]);

yrange = [0.01 0.1 1];
set(gca,'YLim', [0.001 1], 'YTick',yrange,'YTickLabel',sprintfc('%1.0f',yrange*100));
set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;

if saveFig
    savefName = sprintf('ContrastThreshold_vs_Downsampling_%s_noiselevel%1.2f',expName, lateNoiseLevel);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end

%% 5. Plot contrast thresholds vs cone density

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
lowess_span = 0.25;

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
 
% fit current and filtered data
colorsGray = [0.7 0.7 0.7; 0.3 0.3 0.3];

% Fit and plot filtered
y_filtered = weibullFit.ctrthresh(3,:)';
[X3,Y3]    = meshgrid(ones(11,1),log10(x));
Z3 = repmat(log10(y_filtered),1,11);
[meshFit3, gof3] = fit([X3(:) Y3(:)], Z3(:), 'lowess','span',lowess_span);

% Extract single lines for separate ratio's
yFit_filtered = 10.^meshFit3(X3(:,1),Y3(:,1));
R2_filtered   = gof3.rsquare;

% Plot fit and markers 
plot(x, yFit_filtered, 'color', colorsGray(1,:), 'LineWidth', 2); hold all;
scatter(x, y_filtered, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colorsGray(1,:), 'LineWidth',2);

% Fit and plot current
y_current = weibullFit.ctrthresh(2,:)';
[X2,Y2] = meshgrid(ones(11,1),log10(x));
Z2 = repmat(log10(y_current),1,11);
[meshFit2, gof2] = fit([X2(:) Y2(:)], Z2(:), 'lowess','span',lowess_span);

% Extract single lines for separate ratio's
yFit_current = 10.^meshFit2(X2(:,1),Y2(:,1));
R2_current= gof2.rsquare;

% Plot fit and markers 
plot(x, yFit_current, 'color', colorsGray(2,:), 'LineWidth', 2); hold all;
scatter(x, y_current, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor',colorsGray(2,:), 'LineWidth',2);

% Plot absorptions as straight line
y_absorptions = weibullFit.ctrthresh(1,:)';
[f_absorptions, gof1] = fit(x,y_absorptions,'power1');
R2_absorptions = gof1.rsquare; % Coefficient of Determination
plot(x, f_absorptions(x), 'color','k', 'LineWidth', 3); hold all;
scatter(x, y_absorptions, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',2);

% Make plot pretty
box off; axis square;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,...
    'XScale','log', 'YScale', 'log')
xlabel('Cone density (cones/deg^2)','FontSize',20); ylabel('Contrast threshold (%)','FontSize',20)
set(gca, 'XTick',10.^[2:4],'XTickLabel',xlabels_eccen, 'XLim', 10.^[2.5 4.0]);

yrange = [0.001 0.01 0.1 1];
set(gca,'YLim', [0.001 1], 'YTick',yrange,'YTickLabel',sprintfc('%1.1f',yrange*100));
set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
h = findobj(gca,'Type','line');
allDataTypeLabels = {downsamplelbls{:}, 'Filtered by RGC DoG','Cone current','Cone absorptions'};
legend([h(end:-1:1)],allDataTypeLabels, 'Location','bestoutside'); legend boxoff
title(sprintf('Contrast threshold vs cone density'))

% Arrow at 4.5 deg eccentricity
xArrow = [0.345 0.345]; % coordinates normalized to figure size, will change if figure size changes
yArrow = [0.75 0.65];
annotation('textarrow',xArrow,yArrow,'String',['4.5' char(176) ' eccen'], ...
            'FontSize',13,'LineWidth',2)

% Save figure if requested
if saveFig
    savefName = sprintf('ContrastThreshold_vs_Conedensity_%s_noiselevel%1.2f_lowessfit%1.1f', ...
        expName, lateNoiseLevel, lowess_span);
    savefig(fullfile(figurePth,[savefName '.fig']))
    hgexport(gcf,fullfile(figurePth,[savefName '.eps']))
    print(gcf,fullfile(figurePth,[savefName '.png']), '-dpng')
end

%% 6. Get cone density and down sample factor at 4.5 deg 
% i.e, the eccentricity of the psychophysical experiment we compare to

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
predictedContrastThreshold = 10.^meshFit(log10(ratioAtIdx),log10(observedConesAtEccen));

%% 7. Make 3D mesh plot

selectDataTypes = 5:nrDataTypes;
lowess_span = 0.25;

% Fit with meshgrid function (using dummy 2D grid)
[X,Y] = meshgrid(log10(downsampleFactors),log10(x));
Z = log10(weibullFit.ctrthresh(selectDataTypes,:))';
[meshFit, gof] = fit([X(:) Y(:)], Z(:), 'lowess','span',lowess_span);


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
colorsRetina = {'r', 'b', 'g', 'k'}; % nasal, superior, temporal, inferior retina
zlift        = 2*[0.01, 0.01, 0.01, 0.01]; % lift markers a tiny bit for visibility

for jj = 1:4
    scatter3(log10(ratioAtIdx(jj)),log10(observedConesAtEccen(jj)), ...
        log10(predictedContrastThreshold(jj))+zlift(jj), 300, ...
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

% Rotate view: Ratio left / cone density right
set(gca, 'View', [-45.6000   11.2000])
set(gca, 'ydir', 'reverse')

% Save figure again
if saveFig
    fName = sprintf('3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitLowess%1.2f_withDots_view2',lowess_span);
    hgexport(fH4, fullfile(figurePth, [fName '.eps']))
    savefig(fH4, fullfile(figurePth, [fName '.fig']))
    print(fH4, fullfile(figurePth, [fName '.png']), '-dpng')
end

%% 8. Plot asymmetries

% Convert predicted thresholds from log10 fraction to fraction
% Retinal coords: % nasal, superior, temporal, inferior
prediction_retina.rgc.cThresholds.mean = predictedContrastThreshold;
prediction_retina.cones.cThresholds.mean = f_absorptions(observedConesAtEccen); 
prediction_retina.current.cThresholds.mean = 10.^meshFit2([1,2,3,4],log10(observedConesAtEccen))';

% Convert thresholds to sensitivity
prediction_retina.rgc.sensitivity.mean = 1./prediction_retina.rgc.cThresholds.mean;
prediction_retina.cones.sensitivity.mean = 1./prediction_retina.cones.cThresholds.mean;
prediction_retina.current.sensitivity.mean = 1./prediction_retina.current.cThresholds.mean;

%% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
averageConeDensity_stimeccen = mean(observedConesAtEccen);

for ii = 1:4
    errorRatioConeDensity(ii) = 2*abs(diff([observedConesAtEccen(ii),averageConeDensity_stimeccen]));
end

% Nasal, superior, temporal, inferior retina
% upper = - doubling diff in cone density from the mean
% lower = + doubling diff in cone density from the mean

prediction_retina.rgc.cThresholds.error.upper = 10.^meshFit(log10(ratioAtIdx),log10(observedConesAtEccen'-errorRatioConeDensity));
prediction_retina.rgc.cThresholds.error.lower = 10.^meshFit(log10(ratioAtIdx),log10(observedConesAtEccen'+errorRatioConeDensity));

prediction_retina.cones.cThresholds.error.upper = f_absorptions(observedConesAtEccen'-errorRatioConeDensity);
prediction_retina.cones.cThresholds.error.lower = f_absorptions(observedConesAtEccen'+errorRatioConeDensity);

prediction_retina.current.cThresholds.error.upper = 10.^meshFit2([1 1 1 1],log10(observedConesAtEccen'-errorRatioConeDensity));
prediction_retina.current.cThresholds.error.lower = 10.^meshFit2([1 1 1 1],log10(observedConesAtEccen'+errorRatioConeDensity));


% Convert error mRGC thresholds to sensitivity
prediction_retina.rgc.sensitivity.error.upper = 1./prediction_retina.rgc.cThresholds.error.upper;
prediction_retina.rgc.sensitivity.error.lower = 1./prediction_retina.rgc.cThresholds.error.lower;

prediction_retina.cones.sensitivity.error.upper = 1./prediction_retina.cones.cThresholds.error.upper;
prediction_retina.cones.sensitivity.error.lower = 1./prediction_retina.cones.cThresholds.error.lower;

prediction_retina.current.sensitivity.error.upper = 1./prediction_retina.current.cThresholds.error.upper;
prediction_retina.current.sensitivity.error.lower = 1./prediction_retina.current.cThresholds.error.lower;

% Convert retinal coords into visual coords: HVM, UVM (inferior), LVM (superior)
r2VF_wMeanHorz  = @(data) [mean([data(1),data(3)]),data(4),data(2)]; % retina to visual field, where we average nasal/retina
vf2r            = @(data) [data(1),data(4),data(3),data(2)];         % visual field to retina, no averaging
vf2vf_wMeanHorz = @(data) [mean([data(1),data(3)]),data(2),data(4)]; % visual field to visual field where we average nasal/retina

prediction_visualfield.rgc.sensitivity.mean_wHorz = r2VF_wMeanHorz(prediction_retina.rgc.sensitivity.mean');
prediction_visualfield.rgc.sensitivity.error.upper_wHorz = r2VF_wMeanHorz(prediction_retina.rgc.sensitivity.error.upper');
prediction_visualfield.rgc.sensitivity.error.lower_wHorz = r2VF_wMeanHorz(prediction_retina.rgc.sensitivity.error.lower');

prediction_visualfield.cones.sensitivity.mean_wHorz = r2VF_wMeanHorz(prediction_retina.cones.sensitivity.mean');
prediction_visualfield.cones.sensitivity.error.upper_wHorz = r2VF_wMeanHorz(prediction_retina.cones.sensitivity.error.upper');
prediction_visualfield.cones.sensitivity.error.lower_wHorz = r2VF_wMeanHorz(prediction_retina.cones.sensitivity.error.lower');

prediction_visualfield.current.sensitivity.mean_wHorz = r2VF_wMeanHorz(prediction_retina.current.sensitivity.mean');
prediction_visualfield.current.sensitivity.error.upper_wHorz = r2VF_wMeanHorz(prediction_retina.current.sensitivity.error.upper');
prediction_visualfield.current.sensitivity.error.lower_wHorz = r2VF_wMeanHorz(prediction_retina.current.sensitivity.error.lower');

%% Get observed behavior + error
% OBSERVED Left HM, UVM, Right HM, LVM
% (data from baseline experiment Himmelberg, Winawer, Carrasco, 2020, JoV)
observed_visualfield.sensitivity.mean  = [46.4938; 28.9764; 47.7887; 34.3813]; % contrast senstivity (%)
observed_visualfield.sensitivity.error = [2.66468; 1.6445; 1.8450; 2.0505];    % contrast senstivity (%)
observed_retina.sensitivity.mean       = vf2r(observed_visualfield.sensitivity.mean); %  L/R, LVM == superior retina, L/R, UVM == inferior retina

observed_visualfield.sensitivity.mean_wHorz  = vf2vf_wMeanHorz(observed_visualfield.sensitivity.mean);
observed_visualfield.sensitivity.error_wHorz = vf2vf_wMeanHorz(observed_visualfield.sensitivity.error);


%% HVA VMA calc
HVAmean.obs         = hva(observed_retina.sensitivity.mean);
VMAmean.obs         = vma(observed_retina.sensitivity.mean);

HVAmean.predRGC     = hva(prediction_retina.rgc.sensitivity.mean);
VMAmean.predRGC     = vma(prediction_retina.rgc.sensitivity.mean);

HVAmean.predCones   = hva(prediction_retina.cones.sensitivity.mean);
VMAmean.predCones   = vma(prediction_retina.cones.sensitivity.mean);

HVAmean.predCurrent = hva(prediction_retina.current.sensitivity.mean);
VMAmean.predCurrent = vma(prediction_retina.current.sensitivity.mean);

HVAerror.obs       = HVAmean.obs + [-6.90, 6.90]; %  From Himmelberg et al. (2020) 
VMAerror.obs       = VMAmean.obs + [-5.65,5.65];  %  From Himmelberg et al. (2020)

HVAerror.predCones = [hva(prediction_retina.cones.sensitivity.error.lower), hva(prediction_retina.cones.sensitivity.error.upper)];
VMAerror.predCones = [vma(prediction_retina.cones.sensitivity.error.lower), vma(prediction_retina.cones.sensitivity.error.upper)];

HVAerror.predRGC   = [hva(prediction_retina.rgc.sensitivity.error.lower), hva(prediction_retina.rgc.sensitivity.error.upper)];
VMAerror.predRGC   = [vma(prediction_retina.rgc.sensitivity.error.lower), vma(prediction_retina.rgc.sensitivity.error.upper)];

HVAerror.predCurrent = [hva(prediction_retina.current.sensitivity.error.lower), hva(prediction_retina.current.sensitivity.error.upper)];
VMAerror.predCurrent = [vma(prediction_retina.current.sensitivity.error.lower), vma(prediction_retina.current.sensitivity.error.upper)];


combHVA = [HVAmean.predCones, HVAmean.predCurrent, HVAmean.predRGC, HVAmean.obs];
combVMA = [VMAmean.predCones, VMAmean.predCurrent, VMAmean.predRGC, VMAmean.obs];

errorCombHVA = [HVAerror.predCones(1), HVAerror.predCurrent(1), HVAerror.predRGC(1), HVAerror.obs(1); ...
                HVAerror.predCones(2), HVAerror.predCurrent(2), HVAerror.predRGC(2), HVAerror.obs(2)];
            
errorCombVMA = [VMAerror.predCones(1), VMAerror.predCurrent(1), VMAerror.predRGC(1), VMAerror.obs(1); ...
                VMAerror.predCones(2), VMAerror.predCurrent(2), VMAerror.predRGC(2), VMAerror.obs(2)];


fprintf('HVA predicted for cones \t %1.2f%% \n',HVAmean.predCones)
fprintf('VMA predicted for cones \t %1.2f%% \n', VMAmean.predRGC)

fprintf('HVA predicted for current \t %1.2f%% \n',HVAmean.predCurrent)
fprintf('VMA predicted for current \t %1.2f%% \n',VMAmean.predCurrent)

fprintf('HVA predicted for mRGC \t\t %1.2f%% \n', HVAmean.predRGC)
fprintf('VMA predicted for mRGC \t\t %1.2f%% \n', VMAmean.predRGC)

fprintf('HVA observed in behavior \t %1.2f%% \n', HVAmean.obs)
fprintf('VMA observed in behavior \t %1.2f%% \n', VMAmean.obs)


%% Bar plot to compare predictions against behavior

condNames   = {'HVM', 'UVM','LVM'};
condColor   = [63, 121, 204; 0 206 209; 228, 65, 69;150 123 182]/255;
yl          = [0 4];
fH5 = figure(5); set(fH5, 'position',[139,244,1373,543], 'color', 'w'); clf; hold all;

% Plot prediction for cones
subplot(151)
bar(1:3, prediction_visualfield.cones.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(1,:)); hold on
errorbar(1:3,prediction_visualfield.cones.sensitivity.mean_wHorz,...
             (prediction_visualfield.cones.sensitivity.mean_wHorz-prediction_visualfield.cones.sensitivity.error.lower_wHorz), ...
             (prediction_visualfield.cones.sensitivity.error.upper_wHorz-prediction_visualfield.cones.sensitivity.mean_wHorz),...
             '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction Cones');

% Plot prediction for current
subplot(152)
bar(1:3, prediction_visualfield.current.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,prediction_visualfield.current.sensitivity.mean_wHorz,...
             (prediction_visualfield.current.sensitivity.mean_wHorz-prediction_visualfield.current.sensitivity.error.lower_wHorz), ...
             (prediction_visualfield.current.sensitivity.error.upper_wHorz-prediction_visualfield.current.sensitivity.mean_wHorz),...
             '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction Current');


% Plot prediction for mRGCs
subplot(153)
bar(1:3, prediction_visualfield.rgc.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(3,:)); hold on
errorbar(1:3,prediction_visualfield.rgc.sensitivity.mean_wHorz,...
             (prediction_visualfield.rgc.sensitivity.mean_wHorz-prediction_visualfield.rgc.sensitivity.error.lower_wHorz), ...
            (prediction_visualfield.rgc.sensitivity.error.upper_wHorz-prediction_visualfield.rgc.sensitivity.mean_wHorz), ...
            '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction mRGCs');

% Plot prediction for behavior from Himmelberg, Winawer, Carrasco 2020 JoV
subplot(154)
bar(1:3, observed_visualfield.sensitivity.mean_wHorz,'EdgeColor','none','facecolor',condColor(4,:)); hold on
errorbar(1:3,observed_visualfield.sensitivity.mean_wHorz,observed_visualfield.sensitivity.error_wHorz, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^yl, 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Behavior');

% Plot Asymmetries in percent
subplot(155); hold on; cla
x_bar = [0.5, 1, 1.5 2; 3 3.5 4 4.5];
for ii = 1:4
    bar(x_bar(:,ii), [combHVA(ii); combVMA(ii)], 0.2, 'EdgeColor','none','facecolor',condColor(ii,:)); hold on
end
errorbar(x_bar(1,:), combHVA, combHVA-errorCombHVA(1,:), errorCombHVA(2,:)-combHVA, '.','color', 'k', 'LineWidth',2);
errorbar(x_bar(2,:), combVMA, combVMA-errorCombVMA(1,:), errorCombVMA(2,:)-combVMA, '.','color', 'k', 'LineWidth',2);
set(gca,'Xlim',[0,5],'Ylim',[-40 60], 'TickDir', 'out', 'XTick', [1, 2.5], ...
    'XTickLabel', {'HVA', 'VMA'}, 'FontSize', 14, 'YScale', 'linear');
box off;ylabel('Asymmetry (%)');  title('Asymmetry');



if saveFig
    hgexport(fH5, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.eps'))
    savefig(fH5, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.fig'))
    print(fH5, fullfile(figurePth, 'Sensitivity_Model_vs_Behavior_4_5eccen_lateNoiseRGCModel.png'), '-dpng')
end