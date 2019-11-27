function [fH13, fH14, fH15, fH16] = visualizeV1CMFHCP(CMFV1_HCP, meridianDataCMF_RV79, ...
                CMF_HH91, eccDeg, fH12, saveFigures, figureDir)

% ------------ Plot CMF vs eccen for Rovamo & Virsu vs Horton & Hoyt model ------------

titleStr = 'CMF Rovamo Virsu 1979 - in retinal coords';
fH13 = plotMeridiansVsEccen(meridianDataCMF_RV79, eccDeg, titleStr, [], saveFigures);

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
fH14 = plotHVAandVMA(meridianDataCMF_RV79, eccDeg, titleStr, saveFigures);


% ------------ Plot HVA and VMA vs eccen for HCP mean subject -----------
fH15 = figure(fH14);
eccen = mean([0,3.5; 1,6; 3.5,7],2);

% Plot HCP integral data points
plot(eccen(1), mean(CMFV1_HCP.hvaAll035), 'co', 'LineWidth', 4);
plot(eccen(1), mean(CMFV1_HCP.vmaAll035), 'mo', 'LineWidth', 4);

plot(eccen(2), mean(CMFV1_HCP.hvaAll16), 'co', 'LineWidth', 4);
plot(eccen(2), mean(CMFV1_HCP.vmaAll16), 'mo', 'LineWidth', 4);

plot(eccen(3), mean(CMFV1_HCP.hvaAll357), 'co', 'LineWidth', 4);
plot(eccen(3), mean(CMFV1_HCP.vmaAll357), 'mo', 'LineWidth', 4);

% Plot 0 line
plot(0:8, zeros(1,9), 'k')

xlim([0, 40])
ylim([-80, 80])

ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
xlabel('Eccentricity (deg)')
titleStr = 'HVA VMA V1 CMF HCP mean subject - SfN Poster';
title(titleStr);

legend({'HVA Rovamo & Virsu (1979)', ...
    'VMA Rovamo & Virsu (1979)', '', ...
    'HCP HVA integral +/- 15, 0-3.5 deg eccen', ...
    'HCP VMA integral +/- 15, 0-3.5 deg eccen', ...
    'HCP HVA integral +/- 15, 1-6 deg eccen', ...
    'HCP VMA integral +/- 15, 1-6 deg eccen', ...
    'HCP HVA integral +/- 15, 3.5-7 deg eccen', ...
    'HCP VMA integral +/- 15, 3.5-7 deg eccen'}, 'Location', 'SouthEast');
legend boxoff;

set(gca,'FontSize', 14,'TickDir', 'out')

if saveFigures
    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH15, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
    
end


% ------------ Plot HVA and VMA vs eccen for HCP mean subject and RGC and Cones and Watson -----------
fH16 = figure; clf; set(gcf, 'Color', 'w', 'Position', [809   762   751   576]); hold all;
h12 = get(fH12, 'Children');

for ii = [15,14,12:-1:7]
    plot(h12(2).Children(ii).XData,h12(2).Children(ii).YData, ...
        'Color', h12(2).Children(ii).Color, ...
        'LineWidth', h12(2).Children(ii).LineWidth, ...
        'LineStyle', h12(2).Children(ii).LineStyle); hold on
end
xlim([0, 40])
eccen = mean([0,3.5; 1,6; 3.5,7],2);

% Plot HCP integral data points
plot(eccen(1), mean(CMFV1_HCP.hvaAll035), 'co', 'LineWidth', 4);
plot(eccen(1), mean(CMFV1_HCP.vmaAll035), 'mo', 'LineWidth', 4);

plot(eccen(2), mean(CMFV1_HCP.hvaAll16), 'co', 'LineWidth', 4);
plot(eccen(2), mean(CMFV1_HCP.vmaAll16), 'mo', 'LineWidth', 4);

plot(eccen(3), mean(CMFV1_HCP.hvaAll357), 'co', 'LineWidth', 4);
plot(eccen(3), mean(CMFV1_HCP.vmaAll357), 'mo', 'LineWidth', 4);

% Plot 0 line
plot(0:40, zeros(1,41), 'k')

xlim([0, 40])
ylim([-80, 80])

ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
xlabel('Eccentricity (deg)')
titleStr = 'HVA VMA V1 CMF HCP mean subject vs Cones vs mRGC RF';
title(titleStr); set(gca,'FontSize', 14,'TickDir', 'out')


legend({'HVA mRGC RF density Watson 2014 - SfN Poster', 'VMA mRGC RF density Watson 2014 - SfN Poster',...
    'HVA Cones Curcio et al (1990) - rgcDisplacement map', 'VMA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'HVA mRGC RF - rgcDisplacement map', 'VMA mRGC RF - rgcDisplacement map', ...
    'HCP HVA integral +/- 15, 0-3.5 deg eccen', ...
    'HCP VMA integral +/- 15, 0-3.5 deg eccen', ...
    'HCP HVA integral +/- 15, 1-6 deg eccen', ...
    'HCP VMA integral +/- 15, 1-6 deg eccen', ...
    'HCP HVA integral +/- 15, 3.5-7 deg eccen', ...
    'HCP VMA integral +/- 15, 3.5-7 deg eccen'}, 'Location', 'SouthEast');
legend boxoff;

if saveFigures
    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH16, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
end


