% s_plot3DThresholdDensityRatio

% Define data params
ratios      = 1:5; % mRGC 2 Cone ratios
expName     = 'conedensity';
subFolder   = 'average';
fitTypeName =  'powerlaw';
colors      = parula(length(ratios)+1);
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
    
    load(fullfile(dataFolder, sprintf('varThresh_rgcResponse_Cones2RGC%d_absorptionrate_13_conedensity', ratio)), 'varThresh');
    errThresh(ratio,:) = varThresh.*100;
    
    if strcmp(fitTypeName, 'powerlaw')
        yThresh = cell2mat(fit.ctrthresh).*100;
        lm = fitlm(log10(xThresh),log10(yThresh));
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
    plot(xThresh, 10.^Z(ii,:), '-', 'color', colors(ii,:), 'LineWidth',4)
end
for ii = ratios
    errorbar(xThresh,  allData(ii,:), errThresh(ii,:), 'Color', 'k', 'LineStyle','none', 'LineWidth', 1);
    scatter(xThresh, allData(ii,:), 80,'MarkerFaceColor', colors(ii,:), 'MarkerEdgeColor','k', 'LineWidth',1)
end

h = findobj(gca,'Type','line');
legend([h(end:-1:1)],labels, 'Location','Best');
legend boxoff;
set(gca, 'XScale', 'log','YScale','log', 'FontSize', 14, 'LineWidth', 2, 'TickDir', 'out')
ylim([0.5, 20]);
xlabel('Cone array density (cells/deg^2)');
ylabel('Contrast thresholds (%)');
box off; grid on;

if saveFigs
    hgexport(fH1,fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'))
    savefig(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'))
    print(fH1, fullfile(figureFolder, 'ContrastThreshold-vs-Density_withPowerlawFit'), '-dpng')
end
%% Plot 3D data

% make x,y grid
[X,Y] = meshgrid(ratios,xThresh);

c = hot(255);

fH2 = figure(99); set(gcf, 'Position', [394    42   881   756], 'Color', 'w'); clf;
surf(X,Y,allData'); colormap(c(4:end,:));

xlabel('RGC:cone ratio', 'FontSize', 20)
ylabel('Cone array density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20);
title('Data: Effect of RGC filtering on contrast threshold')

set(gca, 'ZLim',[0,16], 'FontSize', 20, 'LineWidth', 2, 'YScale', 'log', 'XScale', 'linear','ZScale', 'log', ...
    'XTick', ratios, 'XTickLabel', sprintfc('%1.2f', 2./ratios), ...
    'YTick', 10.^[2:5], 'YTickLabel', sprintfc('10^%d', [2:5]), ...
    'View',[-152.4000   14.4000], 'TickDir', 'out');
axis square; lighting gouraud

if saveFigs
    hgexport(fH2,fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
    savefig(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'))
    print(fH2, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_Data'), '-dpng')
end

%% Plot fit to data

% upsample ratios
ratioQ     = linspace(1,5,41);
ZUpsampled = interp1(ratios, 10.^Z, ratioQ);
[X,Y] = meshgrid(ratioQ,xThresh);

fH3 = figure(100); set(gcf, 'Position', [  782    44   881   756], 'Color', 'w'); clf;
surf(X,Y,ZUpsampled', 'FaceLighting', 'gouraud')%, 'AmbientStrength', 0.4, 'DiffuseStrength', 0.5, 'SpecularStrength', 0.2, 'SpecularExponent', 40);
colormap(c(20:end,:));
axis square; material shiny;
% lightangle(133,-2)

xlabel('mRGC:cone ratio', 'FontSize', 20)
ylabel('Cone density (cells/deg^2)', 'FontSize', 20)
zlabel('Contrast threshold (%)', 'FontSize', 20)
title('Interpolated fit to data: Effect of RGC filtering on contrast threshold')

set(gca, 'ZLim',[0,10], 'FontSize', 20, 'LineWidth', 2, 'YScale', 'log', 'ZScale', 'log', ...
    'XTick', ratios, 'XTickLabel', sprintfc('%1.2f', 2./ratios), ...
    'View',[-133.2000   25.6000], 'TickDir', 'out')

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
lm_rn = fitlm(log10(xThresh), log10(ZUpsampled(rn,:)));
lm_rs = fitlm(log10(xThresh), log10(ZUpsampled(rs,:)));
lm_rt = fitlm(log10(xThresh), log10(ZUpsampled(rt,:)));
lm_ri = fitlm(log10(xThresh), log10(ZUpsampled(ri,:)));


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
predictedError = [];
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
predictedError(4,2) = cThreshold(observedConesAtEccen(4)+errorRatioConeDensity(4), lm_ri.Coefficients.Estimate(2), lm_ri.Coefficients.Estimate(1));


%% Plot observed/biological variations in cone:mRGC ratios on mesh

figure(fH3); hold all;
colorsRetina = {'r', 'b', 'g', 'k'};
for jj = 1:4
    scatter3(2./ratio_at_idx(jj),observedConesAtEccen(jj),10.^predictedContrastThreshold(jj), 50, 'MarkerFaceColor', colorsRetina{jj}, 'MarkerEdgeColor','w', 'LineWidth',2)
end

% Save figure 3 again with dots
if saveFigs
    hgexport(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots'))
    savefig(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots'))
    print(fH3, fullfile(figureFolder, '3Dmesh_Ratio-vs-Density-vs-Threshold_loglogFitUpsampled_withDots'), '-dpng')
end

%% Make gif of rotating mesh

filename = fullfile(figureFolder,'fullModelMesh.gif');

x = -160:5:-95;
angles = [x, fliplr(x)];

figure(fH3); axis vis3d;
for n = 1:length(angles)
    set(gca, 'View', [angles(n), 10.4000]);
   
    drawnow; pause(0.1);
    
    % Capture the plot as an image
    frame = getframe(fH3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
    end
end



%%

%
% allOptimalRatios = [];
% for r = 1:5
%
%     dataFolder = fullfile(baseFolder, 'data',  'conedensity', 'amplitudes');
%     load(fullfile(dataFolder, sprintf('amplitudesAtStim_ratio1_5_run%d', r)), 'ampsStim');
%
%     allOptimalRatios = cat(4, allOptimalRatios, ampsStim);
% end
% averageOptimalRatio = nanmean(allOptimalRatios,4);
%
% for c = 1:15
%     upsampleAverageOptimalRatio(c,:,:) = interp1(ratios, squeeze(averageOptimalRatio(c,:,:))', ratioQ)';
% end
%
% figure(101); clf;hold all;
% imagesc(squeeze(max(upsampleAverageOptimalRatio))'); colorbar;
% ylabel('ratios')
% xlabel('eccentricities')
%
% [maxRatio,maxRatioIdx] = max(squeeze(max(upsampleAverageOptimalRatio)));
% optimalLine = []
% for ii = 1:length(maxRatioIdx);
%     optimalLine(ii) = xThresh(maxRatioIdx(ii)');
% end
%
% figure(101); hold on;
% plot(maxRatioIdx,1:41, 'r*')
%
% figure(fH3); hold all;
% plot3(ratioQ,optimalLine, ZUpsampled(maxRatioIdx), 'MarkerSize', 500, 'MarkerFaceColor', 'r', 'LineWidth',2)



%% COMPARE TO MODEL TO BEHAVIOR

% Convert thresholds to sensitivity
predContrastSensitivityMEAN    = 1./((10.^predictedContrastThreshold)./100);
predContrastSensitivityMEAN_wHorz_VF = [mean(predContrastSensitivityMEAN([1 3])); predContrastSensitivityMEAN(4); predContrastSensitivityMEAN(2)];
% predContrastSensitivityERROR_upper    = 1./(10.^predictedError(:,2));
predContrastSensitivityERROR_lower    = 1./((10.^predictedError(:,1))./100);

predError = predContrastSensitivityMEAN-predContrastSensitivityERROR_lower';

predContrastSensitivityERROR_wHorz_VF = [mean(predError([1 3])); predError(4); predError(2)];

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
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[1.3, 1.7], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; ylabel('Contrast sensitivity (%)'); title('Model Prediction');

subplot(122)
bar(1:3, obsContrastSensitivityMEAN_wHorz_VF,'EdgeColor','none','facecolor',condColor(2,:)); hold on
errorbar(1:3,obsContrastSensitivityMEAN_wHorz_VF,obsContrastSensitivityERROR_wHorz_VF, '.','color', 'k', 'LineWidth',2);

% format figure
set(gca,'Xlim',[0.2,3.8],'Ylim',10.^[1.3, 1.7], 'TickDir', 'out', 'XTick', [1:3], ...
    'XTickLabel', condNames, 'FontSize', 14, 'YScale', 'log');
box off; title('Behavior');

hva(1./predContrastSensitivityMEAN)
hva(1./obsContrastSensitivityMEAN_VF)

vma(1./predContrastSensitivityMEAN)
vma(1./obsContrastSensitivityMEAN_VF)


if saveFigs
    hgexport(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    savefig(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'))
    print(fH4, fullfile(figureFolder, 'Sensitivity_Model_vs_Behavior_4_5eccen'), '-dpng')
end
