function [fH13, fH14, fH15] = visualizeV1CMFHCP(CMFV1_HCP, meridianDataCMF_RV79, ...
                CMF_HH91, eccDeg, fH10, saveFigures, figureDir)

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
fH14 = plotHVAandVMA(meridianDataCMF_RV79, eccDeg, titleStr, figureDir, saveFigures);




% ------------ Plot HVA and VMA vs eccen for HCP mean subject -----------
fH15 = figure(); clf; set(gcf, 'Color', 'w', 'Position', [ 788   676   754   669]); hold all;

% Boostrap median CMF
nboot = 1000;
bootDataUpperVisualField = bootstrp(nboot, @(x) nanmedian(x), allUprCMF);
bootDataLowerVisualField = bootstrp(nboot, @(x) nanmedian(x), allLowrCMF);
bootDataHorizontalVisualField = bootstrp(nboot, @(x) nanmedian(x), allHorzCMF);
bootDataVerticalVisualField = bootstrp(nboot, @(x) nanmedian(x), allVertCMF);

% Get the median across bootstraps
medianUpperVF = median(bootDataUpperVisualField,1);
medianLowerVF = median(bootDataLowerVisualField,1);
medianHorzVF = median(bootDataHorizontalVisualField,1);
medianVertVF = median(bootDataVerticalVisualField,1);

seUpperVF = std(bootDataUpperVisualField,[],1);
seLowerVF = std(bootDataLowerVisualField,[],1);
seHorzVF = std(bootDataHorizontalVisualField,[],1);
seVertVF = std(bootDataVerticalVisualField,[],1);

% Get eccentricity (to fit data and prepare x-axes for visualization)
eccenHCP = mean([1, 2; ... (1.5 degree)
    2, 3; ... (2.5 degree)
    3, 4; ... (3.5 degree)
    4, 5; ... (4.5 degree)
    5, 6],2); % (5.5 degree)

    % Fit a linear function to data
% Note horizontal VF = nasal+temporal retina (1st/3rd row)
x = (1:0.05:6);
[coeffHorzVF]    = fit(eccenHCP,medianHorzVF','power1');
fitHorzVFHCP     = coeffHorzVF.a .* x.^coeffHorzVF.b;

% Note lower VF = superior retina (2th row)
[coeffLowerVF]    = fit(eccenHCP,medianLowerVF','power1');
fitLowerVFHCP     = coeffLowerVF.a .* x.^coeffLowerVF.b;

% Note upper VF = inferior retina (4th row)
[coeffUpperVF]    = fit(eccenHCP,medianUpperVF','power1');
fitUpperVFHCP     = coeffUpperVF.a .* x.^coeffUpperVF.b;

% Plot HCP integral data points
for ii =  1:length(eccenHCP)
    mdAsym = nanmedian(bootData(ii,:));
    stdAsym = std(bootData(ii,:), 'omitnan');
    
    if mod(ii,2)
        subplot(2,1,1); hold on;
        plot(allEccen(ii), mdAsym, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 12);
        errorbar(allEccen(ii), mdAsym, stdAsym,'Color', 'k', 'LineWidth', lw+1);
    else
        subplot(2,1,2); hold on;
        plot(allEccen(ii), mdAsym, 'MarkerFaceColor', 'k', 'Color', 'k', 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 12);
        errorbar(allEccen(ii), mdAsym, stdAsym, 'Color', 'k', 'LineWidth', lw+1);
    end
end


plot([0 40],[0 0], 'k')
legend boxoff;

axes(fH15.Children(2)); 
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', [-20,130], 'TickDir', 'out', 'FontSize', 14)
title('HVA V1/V2 surface area')

axes(fH15.Children(3)); 
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', [-20,130], 'TickDir', 'out', 'FontSize', 14)
title('VMA V1/V2 surface area')

savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_V1_CMF'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_V1_CMF'), '-dpdf', '-fillpage')


% ------------ Plot HVA and VMA vs eccen for HCP mean subject and RGC and Cones and Watson -----------
% fH16 = figure; clf; set(gcf, 'Color', 'w', 'Position', [809   762   751   576]); hold all;
% h12 = get(fH12, 'Children');
% 
% for ii = [15,14,12:-1:7]
%     plot(h12(2).Children(ii).XData,h12(2).Children(ii).YData, ...
%         'Color', h12(2).Children(ii).Color, ...
%         'LineWidth', h12(2).Children(ii).LineWidth, ...
%         'LineStyle', h12(2).Children(ii).LineStyle); hold on
% end
% xlim([0, 40])
% eccen = mean([0,3.5; 1,6; 3.5,7],2);
% 
% % Plot HCP integral data points
% plot(eccen(1), mean(CMFV1_HCP.hvaAll035), 'co', 'LineWidth', 4);
% plot(eccen(1), mean(CMFV1_HCP.vmaAll035), 'mo', 'LineWidth', 4);
% 
% plot(eccen(2), mean(CMFV1_HCP.hvaAll16), 'co', 'LineWidth', 4);
% plot(eccen(2), mean(CMFV1_HCP.vmaAll16), 'mo', 'LineWidth', 4);
% 
% plot(eccen(3), mean(CMFV1_HCP.hvaAll357), 'co', 'LineWidth', 4);
% plot(eccen(3), mean(CMFV1_HCP.vmaAll357), 'mo', 'LineWidth', 4);
% 
% % Plot 0 line
% plot(0:40, zeros(1,41), 'k')
% 
% xlim([0, 40])
% ylim([-80, 80])
% 
% ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
% xlabel('Eccentricity (deg)')
% titleStr = 'HVA VMA V1 CMF HCP mean subject vs Cones vs mRGC RF';
% title(titleStr); set(gca,'FontSize', 14,'TickDir', 'out')
% 
% 
% legend({'HVA mRGC RF density Watson 2014 - SfN Poster', 'VMA mRGC RF density Watson 2014 - SfN Poster',...
%     'HVA Cones Curcio et al (1990) - rgcDisplacement map', 'VMA Cones Curcio et al (1990) - rgcDisplacement map', ...
%     'HVA mRGC RF - rgcDisplacement map', 'VMA mRGC RF - rgcDisplacement map', ...
%     'HCP HVA integral +/- 15, 0-3.5 deg eccen', ...
%     'HCP VMA integral +/- 15, 0-3.5 deg eccen', ...
%     'HCP HVA integral +/- 15, 1-6 deg eccen', ...
%     'HCP VMA integral +/- 15, 1-6 deg eccen', ...
%     'HCP HVA integral +/- 15, 3.5-7 deg eccen', ...
%     'HCP VMA integral +/- 15, 3.5-7 deg eccen'}, 'Location', 'SouthEast');
% legend boxoff;
% 
% if saveFigures
%     % Save matlab fig and pdf
%     figName = strrep(titleStr,' ','_');
%     savefig(fH16, fullfile(figureDir, figName))
%     print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
% end
% 

