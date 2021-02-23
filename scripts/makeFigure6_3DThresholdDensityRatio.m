function makeFigure6_3DThresholdDensityRatio()
% Function to make Figure 6 of the manuscript:
%   Radial asymmetries around the visual field: From retina to cortex to 
%   behavior. By Kupers, Benson, Carrasco, Winawer.
%    JOURNAL. DOI.

% This function requires you to have psychometric functions of the 
% simulated data. You can run this with:
% 
% for ratio = 1:5
%     plotPsychometricFunctionsRGCModel(baseFolder, 'conedensity', 'average',ratio, 'plotAvg',true, 'meanPoissonPaddingFlag', true);
% end

%% 0. Define params and folders
nrSimConeDensities = 13;
c2rgc       = 1:5; % Cone 2 RGC ratios
rgc2c       = 2./c2rgc; % RGC 2 Cone ratios
expName     = 'conedensity';
meanPoissonPaddingFlag = true;
colors      = parula(length(c2rgc)+1);
labels      = sprintfc('RGC:cone = %1.1f:1.0', rgc2c);

% Change folder names if using mean Poisson padded cone data
if meanPoissonPaddingFlag
    extraSubFolder = 'withPaddingBeforeConvolution';
    subFolder   = 'average_meanPoissonPadded';
else
    extraSubFolder = 'noPaddingBeforeConvolution';
    subFolder   = 'average';
end

% Folders
saveFigs     = true;
baseFolder   = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
dataFolder   = fullfile(baseFolder,'data',expName,'thresholds',extraSubFolder);
figureFolder = fullfile(baseFolder,'figures','surface3D', expName, subFolder,extraSubFolder);
if (saveFigs) && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

%% 1. Load thresholds for each cone2mrgc ratio and fit function to data

% preallocate space
allData   = NaN(length(c2rgc), nrSimConeDensities);
errThresh = allData;
fct       = allData;
R2        = NaN(length(c2rgc),1);

% Loop over ratios
for r = c2rgc
    
    % Load simulated contrast thresholds
    load(fullfile(dataFolder, sprintf('cThresholds_ratio%d_%s', r, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');
    allData(r,:) = [reshape(cell2mat(fit.ctrthresh),[],length(fit.ctrthresh))].*100;
    
    % Load variance in simulated contrast thresholds computed by
    % bootstrapping simulation iterations starting with different rng seeds
    load(fullfile(dataFolder, sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptionrate_13_conedensity', r)), 'varThresh');
    errThresh(r,:) = varThresh.*100;
    
    % Fit a linear function in log-log space (so a powerlaw) for each simulated mrgc2cone ratio
    contrastThresh = cell2mat(fit.ctrthresh).*100; 
    coneDensities  = xThresh;
    lm = fitlm(log10(coneDensities),log10(contrastThresh));
    
    % Get fitted contrast thresholds (fct)
    fct(r,:) = lm.Coefficients.Estimate(2).*log10(xThresh) + lm.Coefficients.Estimate(1);
    
    % Get R2 of fits
    R2(r,:) = lm.Rsquared.ordinary;
end

% clean up 
clear fit xThresh fitToPlot dataToPlot;

%% 2. Upsample ratios and corresponding thresholds

dt = 41; % amount of upsampled values
c2rgc_upsampled  = linspace(c2rgc(1),c2rgc(end),dt);
rgc2c_upsampled  = 2./c2rgc_upsampled;
% fct_upsampled     = interp1(c2rgc, 10.^fct, ratios_upsampled);
% [X,Y]             = meshgrid(ratios_upsampled,coneDensities);

% Resample to log space, to get rectangles in mesh
coneDensities_resampled = 10.^linspace(log10(coneDensities(1)),log10(coneDensities(end)),13);

[X1,Y1] = meshgrid(c2rgc,coneDensities);
[X2,Y2] = meshgrid(c2rgc_upsampled,coneDensities_resampled);

fct_upsampled2 = griddata(X1,Y1, 10.^fct', X2,Y2, 'linear');

%% ---------------------------------
% ---------- Visualize ----------- %
% ----------------------------------

% Plot threshold fits vs density data for each ratio
fH1 = figure(1); set(gcf, 'Color', 'w', 'Position', [79 39 1522 759]); clf; hold all;

for ii = c2rgc
    subplot(2,3,ii); hold on;
    scatter(coneDensities,  allData(ii,:), 40, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',1)
    errorbar(coneDensities, allData(ii,:), errThresh(ii,:), 'Color', 'k', 'LineStyle','none', 'LineWidth', 1);
    plot(coneDensities, 10.^fct(ii,:), 'r-', 'LineWidth',4)
    
    set(gca, 'XScale', 'log','YScale','log', 'FontSize', 14, 'LineWidth', 2, 'TickDir', 'out')
    ylim([0.5, 20]);
    title(sprintf('%1.1f mRGCs vs 1 cone', 2/ii));
    xlabel('Cone array density (cells/deg^2)');
    ylabel('Contrast thresholds (%)');
    box off; grid on;
end

% Plot lines in additional subplot
subplot(2,3,6); hold all;
for ii = c2rgc
    plot(coneDensities, 10.^fct(ii,:), '-', 'LineWidth',4, 'Color', colors(ii,:))
end

h = findobj(gca,'Type','line');
legend([h(end:-1:1)],labels, 'Location','Best');
legend boxoff;
set(gca, 'XScale', 'log','YScale','log', 'FontSize', 14, 'LineWidth', 2, 'TickDir', 'out')
ylim([0.5, 20]);
title('All 5 mRGC : cones ratios')
xlabel('Cone array density (cells/deg^2)');
ylabel('Contrast thresholds (%)');
box off; grid on;

if saveFigs
    hgexport(fH1,fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'))
    savefig(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'))
    print(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'), '-dpng')
end

%% Plot 3D data without interpolation (RGC2Cone ratio vs Cone Density vs Threshold)

% % Make x,y grid
% [X,Y] = meshgrid(ratios,coneDensities);
%
% % Define color map
% c = gray(255);
%
% fH2 = figure(2); set(gcf, 'Position', [394    42   881   756], 'Color', 'w'); clf;
% surf(X,Y,allData'); colormap(c(4:end,:));
%
% xlabel('RGC:cone ratio', 'FontSize', 20)
% ylabel('Cone array density (cells/deg^2)', 'FontSize', 20)
% zlabel('Contrast threshold (%)', 'FontSize', 20);
% title('Data: Effect of RGC filtering on contrast threshold')
%
% set(gca, 'ZLim',[0,16], 'FontSize', 20, 'LineWidth', 2, 'YScale', 'log', 'XScale', 'linear','ZScale', 'log', ...
%     'XTick', ratios, 'XTickLabel', sprintfc('%1.2f', 2./ratios), ...
%     'YTick', 10.^[2:5], 'YTickLabel', sprintfc('10^%d', [2:5]), ...
%     'View',[-152.4000   14.4000], 'TickDir', 'out');
% axis square; lighting gouraud
%
% if saveFigs
%     hgexport(fH2,fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
%     savefig(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
%     print(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'), '-dpng')
% end

%% Plot 3D mesh with INTERPOLATED data

fH3 = figure(3); clf; set(gcf, 'Position', [782 44 881 756], 'Color', 'w'); clf;
surf(X2,Y2,fct_upsampled2, 'FaceLighting', 'gouraud', 'FaceColor',[1 1 1]);% 'DiffuseStrength', 0.5, 'SpecularStrength', 0.2, 'SpecularExponent', 40);
% camlight(-180,0,'local')
% h = light;
% lightangle(h,90,100)

% Add labels
xlabel('mRGC:cone ratio', 'FontSize', 20)
ylabel('Cone density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Interpolated fit to data: Effect of RGC filtering on contrast threshold')

% Make plot pretty
set(gca, 'ZLim',[0,10], 'FontSize', 20, 'LineWidth', 2, ...
    'XScale', 'linear', 'YScale', 'log', 'ZScale', 'log', ...
    'XTick', c2rgc, 'XTickLabel', sprintfc('%1.2f', rgc2c), ...
    'TickDir', 'out','View',[-134.0000   11.2000]);
grid on;
set(gca, 'GridAlpha', .1, 'ZMinorGrid', 'on', 'YMinorGrid', 'off', 'XMinorGrid', 'off')
set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1])
axis square; material shiny;

% Save figure if  requested
if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view1'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view1'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view1'), '-dpng')
end

% Flip view and save again
set(gca, 'View', [-219.6000   12.0000])

if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view2'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view2'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_view2'), '-dpng')
end


%%  Predict contrast thresholds for given mRGC density

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

% Compute RGC:cone ratio
rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2015.eccDeg==eccToCompute); % index
ratioAtIdx   = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Find the ratio of interest in model grid for nasal and inferior retina
[err_rn,rn] = min(abs(rgc2c_upsampled-ratioAtIdx(1))); % nasal retina
[err_rs,rs] = min(abs(rgc2c_upsampled-ratioAtIdx(2))); % superior retina
[err_rt,rt] = min(abs(rgc2c_upsampled-ratioAtIdx(3))); % temporal retina
[err_ri,ri] = min(abs(rgc2c_upsampled-ratioAtIdx(4))); % inferior retina

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratioAtIdx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% Fit new line to upsampled contrast threshold data for all meridians
fct_upsampled2 = fct_upsampled2';
lm_rn = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rn,:))); % nasal
lm_rs = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rs,:))); % superior
lm_rt = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(rt,:))); % temporal
lm_ri = fitlm(log10(coneDensities_resampled), log10(fct_upsampled2(ri,:))); % inferior

% Predicted threshold
cThreshold = @(x, a_coeff, b_intcpt) (a_coeff* log10(x)) + b_intcpt;

% Predicted cone density
% predictedDensity = @(y, a_coeff, b_intcpt) (10.^((y-b_intcpt)./a_coeff));

% Get contrast levels from model at 4.5 deg eccen:
% Nasal retina
predictedContrastThreshold(1) = cThreshold(observedConesAtEccen(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedContrastThreshold(2) = cThreshold(observedConesAtEccen(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedContrastThreshold(3) = cThreshold(observedConesAtEccen(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedContrastThreshold(4) = cThreshold(observedConesAtEccen(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));

% Define error margins in terms of cone density, i.e.:
% Double (upper bound) or half (lower bound) the difference in cone density from the mean, for each meridian
averageConeDensity_stimeccen = mean(observedConesAtEccen);

for ii = 1:4
    errorRatioConeDensity(ii) = 2*diff([observedConesAtEccen(ii),averageConeDensity_stimeccen]);
end

% Get predicted thresholds for those error margins, using their model fits
predictedError = [];

% Nasal retina
predictedError(1,1) = cThreshold(observedConesAtEccen(1)-errorRatioConeDensity(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
predictedError(1,2) = cThreshold(observedConesAtEccen(1)+errorRatioConeDensity(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedError(2,1) = cThreshold(observedConesAtEccen(2)-errorRatioConeDensity(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
predictedError(2,2) = cThreshold(observedConesAtEccen(2)+errorRatioConeDensity(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedError(3,1) = cThreshold(observedConesAtEccen(3)-errorRatioConeDensity(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
predictedError(3,2) = cThreshold(observedConesAtEccen(3)+errorRatioConeDensity(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedError(4,1) = cThreshold(observedConesAtEccen(4)-errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));
predictedError(4,2) = cThreshold(observedConesAtEccen(4)+errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1)); %#ok<NASGU>

%% Plot observed/biological variations in cone:mRGC ratios on mesh
figure(fH3); hold all;
colorsRetina = {'r', 'b', 'g', 'k'};
zlift        = [0.01, 0.03, 0.01, 0.1]; % lift markers a tiny bit for visibility
for jj = 1:4
    scatter3(2./ratioAtIdx(jj),observedConesAtEccen(jj),10.^predictedContrastThreshold(jj)+zlift(ii), 300, 'MarkerFaceColor', colorsRetina{jj}, 'MarkerEdgeColor','k', 'LineWidth',0.1, 'Marker', 'p')
end

% Save figure 3 again with dots
if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots_view2'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots_view2'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots_view2'), '-dpng')
end

% Ratio left / cone density right
set(gca, 'View', [-134.0000   11.2000])
set(gca, 'xdir', 'reverse')

if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_withDots_view1'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_withDots_view1'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled_withDots_view1'), '-dpng')
end

%% Make gif of rotating mesh
%
% filename = fullfile(figureFolder,'fullModelMesh.gif');
%
% x = -160:5:-95;
% angles = [x, fliplr(x)];
%
% axis manual; axis vis3d;
% for n = 1:length(angles)
%     set(gca, 'View', [angles(n), 10.4000]);
%
%     drawnow; pause(0.1);
%
%     % Capture the plot as an image
%     frame = getframe(fH3);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if n == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
%     end
% end


