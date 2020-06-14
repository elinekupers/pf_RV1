% s_plot3DThresholdDensityRatio

% Define data params
ratios      = 1:5; % mRGC 2 Cone ratios
expName     = 'conedensity';
subFolder   = 'run1';
fitTypeName =  'linear';
colors      = jet(length(ratios));
labels      = sprintfc('RGC:cone = %1.1f:1.0', 2./ratios);

% Folders
saveFigs     = true;
baseFolder   = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
dataFolder   = fullfile(baseFolder,'data',expName,'thresholds');
figureFolder = fullfile(baseFolder,'figures','surface3D', expName, subFolder);
if (saveFigs) && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

allData = [];
for ratio = ratios
    load(fullfile(dataFolder, sprintf('cThresholds_ratio%d_%s', ratio, subFolder)), 'expName','expParams', 'dataToPlot', 'fitToPlot','fit', 'xThresh');
    allData(ratio,:) = [reshape(cell2mat(fit.ctrthresh),[],length(fit.ctrthresh))].*100;
    
    if strcmp(fitTypeName, 'linear')
        yThresh = cell2mat(fit.ctrthresh).*100;
        lm = fitlm(log10(xThresh),yThresh);
        R2(ratio,:) = lm.Rsquared.ordinary;
        
        Z(ratio,:) = lm.Coefficients.Estimate(2).*log10(xThresh) + lm.Coefficients.Estimate(1);
        
        
    elseif strcmp(fitTypeName, 'poly2')
        yThresh = cell2mat(fit.ctrthresh)'.*100;
        [fitResult, err]  = polyfit(log10(xThresh),yThresh,2);
        Z(ratio,:)  = polyval(fitResult,log10( xThresh), err);
        R2 = 1 - (err.normr/norm(yThresh - mean(yThresh)))^2;
    end
    
end

%% Plot 2D data

fH1 = figure(98); set(gcf, 'Color', 'w'); clf; hold all;
for ii = ratios
    scatter(xThresh, allData(ii,:), 50,'MarkerFaceColor', colors(ii,:), 'MarkerEdgeColor','k', 'LineWidth',2)
    plot(xThresh, Z(ii,:), '-', 'color', colors(ii,:), 'LineWidth',2)
end

h = findobj(gca,'Type','line');
legend([h(end:-1:1)],labels, 'Location','Best');
legend boxoff;
set(gca, 'XScale', 'log', 'FontSize', 14, 'LineWidth', 2, 'TickDir', 'out')
xlabel('Cone array density (cells/deg^2)');
ylabel('Contrast thresholds (%)');
box off; grid on;

if saveFigs
    hgexport(fH1,fullfile(figureFolder, 'ContrastThreshold-vs-Density_withLinearFit'))
    savefig(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withLinearFit'))
    print(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withLinearFit'), '-dpng')
end
%% Plot 3D data

% make x,y grid
[X,Y] = meshgrid(ratios,xThresh);

fH2 = figure(99); set(gcf, 'Position', [394    42   881   756], 'Color', 'w'); clf; 
surf(X,Y,allData')

xlabel('RGC:cone ratio', 'FontSize', 20)
ylabel('Cone array density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Data: Effect of RGC filtering on contrast threshold')

set(gca, 'ZLim',[0,16], 'FontSize', 20, 'LineWidth', 2, 'YScale', 'log', ...
         'XTick', ratios, 'XTickLabel', sprintfc('%1.2f', 2./ratios), ...
         'View',[-111.6000   22.4000], 'TickDir', 'out')

if saveFigs
    hgexport(fH2,fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
    savefig(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
    print(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'), '-dpng')
end

%% Plot fit to data

% upsample ratios
ratioQ = linspace(1,5,41);
ZUpsampled = interp1(ratios, Z, ratioQ);
[X,Y] = meshgrid(ratioQ,xThresh);

fH3 = figure(100); set(gcf, 'Position', [394    42   881   756], 'Color', 'w'); clf; 
surf(X,Y,ZUpsampled')

xlabel('RGC:cone ratio', 'FontSize', 20)
ylabel('Cone array density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Interpolated fit to data: Effect of RGC filtering on contrast threshold')

set(gca, 'ZLim',[0,10], 'FontSize', 20, 'LineWidth', 2, 'YScale', 'log', ...
         'XTick', ratios, 'XTickLabel', sprintfc('%1.2f', 2./ratios), ...
         'View',[-111.6000   22.4000], 'TickDir', 'out')

if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_linearFitUpsampled'), '-dpng')
end


%%  Predict contrast thresholds for given mRGC density

% mRGC data for different meridia. Order=nasal, superior, temporal,inferior.
watson2015 = load(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2');
watson2015.eccDeg = (0:0.05:60);

% Curcio et al. 1990 (left eye, retina coords, same order as mRGC)
curcio1990 = load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'conesCurcioIsetbio', 'eccDeg','angDeg');
[~,angIdx] = intersect(curcio1990.angDeg,[0,90,180,270]);
coneDensityDeg2PerMeridian= curcio1990.conesCurcioIsetbio(angIdx,:);

% Compute RGC:cone ratio
rgc2coneRatio = watson2015.mRGCRFDensityPerDeg2./coneDensityDeg2PerMeridian;

% Get rgc:cone ratio at chosen eccentricity
eccToCompute = 4.5; % deg
idxEccen     = find(watson2015.eccDeg==eccToCompute); % index
ratio_at_idx = rgc2coneRatio(:,idxEccen); % mRGC:cone ratio at index

% Find the ratio of interest in model grid for nasal and inferior retina
upsampledRatios = 2./ratioQ; % go from cone:rgc to 2rgc:cone 
[~,rn] = min(abs(upsampledRatios-ratio_at_idx(1))); % nasal retina
[~,rs] = min(abs(upsampledRatios-ratio_at_idx(2))); % superior retina
[~,rt] = min(abs(upsampledRatios-ratio_at_idx(3))); % temporal retina
[~,ri] = min(abs(upsampledRatios-ratio_at_idx(4))); % inferior retina

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratio_at_idx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute))

% Fit new line to upsampled contrast threshold data for nasal and inferior mRGC:cone ratio
lm_rn = fitlm(log10(xThresh), ZUpsampled(rn,:));
lm_rs = fitlm(log10(xThresh), ZUpsampled(rs,:));
lm_rt = fitlm(log10(xThresh), ZUpsampled(rt,:));
lm_ri = fitlm(log10(xThresh), ZUpsampled(ri,:));

% Predicted threshold
cThreshold = @(x, a_coeff, b_intcpt) (a_coeff* log10(x)) + b_intcpt;

% Predicted cone density
predictedDensity = @(y, a_coeff, b_intcpt) (10.^((y-b_intcpt)./a_coeff));

% Get contrast levels from model
predictedContrastThreshold.nasalRetina    = cThreshold(observedConesAtEccen(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
predictedContrastThreshold.superiorRetina = cThreshold(observedConesAtEccen(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
predictedContrastThreshold.temporalRetina = cThreshold(observedConesAtEccen(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
predictedContrastThreshold.inferiorRetina = cThreshold(observedConesAtEccen(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));


