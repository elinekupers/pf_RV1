function [fH13, fH14, fH15] = visualizeV1CMFHCP(V1CMF, meridianDataCMF_RV79, ...
                CMF_HH91, eccDeg, saveFigures, figureDir)

% ------------ Plot CMF vs eccen for Rovamo & Virsu vs Horton & Hoyt model ------------

titleStr = 'CMF Rovamo Virsu 1979 - in retinal coords';
fH13 = plotMeridiansVsEccen(meridianDataCMF_RV79, eccDeg, titleStr, [], figureDir, saveFigures);

figure(fH13); hold on;
plot(eccDeg, CMF_HH91, 'LineWidth', 3, 'Color', [.7 .7 .7])
legend({'nasal meridian on retina',...
    'superior meridian on retina', ...
    'temporal meridian on retina', ...
    'inferior meridian on retina', ...
    'Hoyton & Hort (1991)'}, 'Location', 'Best');

titleStr = 'CMF Rovamo Virsu 1979, Hoyton Hort 1991 - in retinal coords';
title(titleStr);
figName = strrep(titleStr,' ','_');
savefig(fH13, fullfile(figureDir, figName))
print(fullfile(figureDir, figName), '-dpdf', '-fillpage')


% ------------ Plot HVA and VMA vs eccen for Rovamo & Virsu model ------------
titleStr = 'HVA VMA CMF Rovamo Virsu 1979 - in retinal coords';
fH14 = plotHVAandVMA(meridianDataCMF_RV79, [], eccDeg, [], titleStr, figureDir, saveFigures);

% ------------ Plot HVA and VMA vs eccen for HCP mean subject -----------
fH15 = figure(); clf; set(gcf, 'Color', 'w', 'Position', [418, 269, 1905, 872]); hold all;

% Plot HCP integral data points
subplot(1,2,1); hold on;
for ii = 1:length(V1CMF.mdHVA)
    plot(V1CMF.polyfit_eccDeg(ii), V1CMF.mdHVA(ii), 'MarkerFaceColor', 'k', 'Color', 'k', 'Marker', 'o', 'LineWidth', 2, 'MarkerSize', 12);
    errorbar(V1CMF.polyfit_eccDeg(ii), V1CMF.mdHVA(ii), V1CMF.stdHVA(ii),'Color', 'k', 'LineWidth', 2);
end
plot(V1CMF.polyval_eccDeg,V1CMF.lineFitHVA, 'k:', 'lineWidth',2)
plot([0 40],[0 0], 'k')
grid on;
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 8], 'ylim', [-20,130],'TickDir', 'out', 'FontSize', 14)
title('HVA V1/V2 surface area')

subplot(1,2,2); hold on;
for ii = 1:length(V1CMF.mdVMA)
    plot(V1CMF.polyfit_eccDeg(ii), V1CMF.mdVMA(ii), 'MarkerFaceColor', 'k', 'Color', 'k', 'Marker', 'o', 'LineWidth', 2, 'MarkerSize', 12);
    errorbar(V1CMF.polyfit_eccDeg(ii), V1CMF.mdVMA(ii), V1CMF.stdVMA(ii), 'Color', 'k', 'LineWidth', 2);
end
plot(V1CMF.polyval_eccDeg,V1CMF.lineFitVMA, 'k:', 'lineWidth',2)
plot([0 40],[0 0], 'k')
grid on;
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 8], 'ylim', [-20,130], 'TickDir', 'out', 'FontSize', 14)
title('VMA V1/V2 surface area')

savefig(fullfile(figureDir, 'HVA_VMA_V1_CMF'))
print(fullfile(figureDir, 'HVA_VMA_V1_CMF'), '-dpdf', '-fillpage')



