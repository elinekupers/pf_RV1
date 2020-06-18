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
[err_rn,rn] = min(abs(upsampledRatios-ratio_at_idx(1))); % nasal retina
[err_rs,rs] = min(abs(upsampledRatios-ratio_at_idx(2))); % superior retina
[err_rt,rt] = min(abs(upsampledRatios-ratio_at_idx(3))); % temporal retina
[err_ri,ri] = min(abs(upsampledRatios-ratio_at_idx(4))); % inferior retina

% Get cone density at chosen eccentricity for each meridian
observedConesAtEccen = watson2015.mRGCRFDensityPerDeg2(:,idxEccen)./ratio_at_idx;

% Check: should be equal to curcio data
isequal(observedConesAtEccen,coneDensityDeg2PerMeridian(:,curcio1990.eccDeg==eccToCompute));

% Fit new line to upsampled contrast threshold data for all meridians
lm_rn = fitlm(log10(xThresh), ZUpsampled(rn,:));
lm_rs = fitlm(log10(xThresh), ZUpsampled(rs,:));
lm_rt = fitlm(log10(xThresh), ZUpsampled(rt,:));
lm_ri = fitlm(log10(xThresh), ZUpsampled(ri,:));


% Predicted threshold
cThreshold = @(x, a_coeff, b_intcpt) (a_coeff* log10(x)) + b_intcpt;

% Predicted cone density
predictedDensity = @(y, a_coeff, b_intcpt) (10.^((y-b_intcpt)./a_coeff));

% Get contrast levels from model at 4.5 deg eccen
% Nasal retina
predictedContrastThreshold(1) = cThreshold(observedConesAtEccen(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedContrastThreshold(2) = cThreshold(observedConesAtEccen(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedContrastThreshold(3) = cThreshold(observedConesAtEccen(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedContrastThreshold(4) = cThreshold(observedConesAtEccen(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));


% Double the difference in cone density from the mean, for each meridian
errorRatioConeDensity = 0.5*abs(observedConesAtEccen-mean(observedConesAtEccen));

% Nasal retina
predictedError(1) = cThreshold(errorRatioConeDensity(1), lm_rn.Coefficients.Estimate(2), lm_rn.Coefficients.Estimate(1));
% Superior retina
predictedError(2) = cThreshold(errorRatioConeDensity(2), lm_rs.Coefficients.Estimate(2), lm_rs.Coefficients.Estimate(1));
% Temporal retina
predictedError(3) = cThreshold(errorRatioConeDensity(3), lm_rt.Coefficients.Estimate(2), lm_rt.Coefficients.Estimate(1));
% Inferior retina
predictedError(4) = cThreshold(errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));




%% COMPARE TO MODEL TO BEHAVIOR

% Convert thresholds to sensitivity
predContrastSensitivityMEAN    = 1./(predictedContrastThreshold./100);
predContrastSensitivityMEAN_wHorz_VF = [mean(predContrastSensitivityMEAN([1 3])); predContrastSensitivityMEAN(4); predContrastSensitivityMEAN(2)];
predContrastSensitivityERROR    = 1./(predictedError./100);
predContrastSensitivityERROR_wHorz_VF = [mean(predContrastSensitivityERROR([1 3])); predContrastSensitivityERROR(4); predContrastSensitivityERROR(2)];



obsContrastSensitivityMEAN_VF = [46.4938; 47.7887; 28.9764; 34.2813];
obsContrastSensitivityERROR_VF = [2.66468; 1.8450; 1.6445; 2.0505];
obsContrastSensitivityMEAN_wHorz_VF = [mean(obsContrastSensitivityMEAN_VF([1 2])), obsContrastSensitivityMEAN_VF(3), obsContrastSensitivityMEAN_VF(4)];
obsContrastSensitivityERROR_wHorz_VF = [mean(obsContrastSensitivityERROR_VF([1 2])), obsContrastSensitivityERROR_VF(3), obsContrastSensitivityERROR_VF(4)];

% Bar plot to compare against behavior

condNames   = {'HVM', 'UVM','LVM'};
condColor   = [63, 121, 204; 228, 65, 69]/255;

fH4 = figure(4); set(fH4, 'position',[ 824   348   688   439], 'color', 'w'); clf; hold all;

subplot(121)
bar(1:3, predContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(1,:)); hold on
errorbar(1:3,predContrastSensitivityMEAN_wHorz_VF,predContrastSensitivityERROR_wHorz_VF,'.','color', 'k', 'LineWidth',2);

% format figure
set(gca,'xlim',[0.2,3.8],'ylim',[0, 70], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14);
box off; ylabel('Contrast sensitivity'); title('Model'); 

subplot(122)
bar(1:3, obsContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,obsContrastSensitivityMEAN_wHorz_VF,obsContrastSensitivityERROR_wHorz_VF, '.','color', 'k', 'LineWidth',2);

% format figure
set(gca,'xlim',[0.2,3.8],'ylim',[0, 70], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14);
box off; title('Behavior'); 


if saveFigs
    hgexport(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    savefig(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    print(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'), '-dpng')
end
