function [fH13, fH14, fH15, fH16] = visualizeV1CMFHCP(CMFV1_HCP, meridianDataCMF_RV79, ...
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
fH15 = figure(); clf; hold all;
eccen = mean([1, 6; ... 1,2 -> 3.5 deg    (5 degree)
              0, 3.5;... 3,4 -> 1.75 deg  (3.5 degree)
              3.5, 7; ... 5,6 -> 5.25 deg (3.5 degree)
              0.5, 1; ... 7,8 -> 0.75 deg (0.5 degree) DATA ARE TOO NOISY
              1, 2; ... 9,10 -> 1.5 deg   (1 degree)   DATA ARE TOO NOISY
              2, 4; ... 11,12 -> 3 deg    (2 degree)
              4, 6; ... 13,14 -> 5 deg    (2 degree)
              4,8; ...  15,16 -> 6 deg    (4 degree)
              6,8],2); % 17,18 -> 7 deg   (2 degree)
fn = fieldnames(CMFV1_HCP);

allEccen = NaN(1,length(eccen)*2);
allEccen(1:2:end) = eccen;
allEccen(2:2:end) = eccen;

 

colors = [255, 165, 0]./255; lw = 2;
% Plot HCP integral data points
for ii =  1:length(allEccen)
    if any(ii==[1:6,11:18])
        if mod(ii,2)
            plot(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), 'MarkerFaceColor', colors, 'MarkerEdgeColor', colors, 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 15);
            errorbar(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), std(CMFV1_HCP.(fn{ii}))/sqrt(181),'Color', 'k', 'LineWidth', lw+1);
        else
            plot(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), 'MarkerFaceColor', [1 1 1], 'Color', colors, 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 15);
            errorbar(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), std(CMFV1_HCP.(fn{ii}))/sqrt(181),'Color', 'k', 'LineWidth', lw+1);
        end
    end         
end
legend({'HVA','', 'VMA'}, 'Location', 'NorthWest'); 

plot([0 40],[0 0], 'k')
legend boxoff;

ylabel('more vertical/upper VF <- Asymmetry (%) -> more horizontal/lower VF')
xlabel('Eccentricity (deg)')
title(titleStr)
set(gca', 'xlim', [0 8], 'ylim', [-80,80], 'TickDir', 'out', 'FontSize', 14)


figure(fH10); hold on;
for ii =  1:length(allEccen)
    if any(ii==[1:6,11:18])
        if mod(ii,2)
            plot(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), 'MarkerFaceColor', colors, 'MarkerEdgeColor', colors, 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 15);
            errorbar(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), std(CMFV1_HCP.(fn{ii}))/sqrt(181),'Color', 'k', 'LineWidth', lw+1);
        else
            plot(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), 'MarkerFaceColor', [1 1 1], 'Color', colors, 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 15);
            errorbar(allEccen(ii), mean(CMFV1_HCP.(fn{ii})), std(CMFV1_HCP.(fn{ii}))/sqrt(181),'Color', 'k', 'LineWidth', lw+1);
        end
    end         
end

obj = findobj(gca,'Type','Line');
legend(obj([21,20,18,17,16,15,12,11]), {'HVA mRGCf density Watson 2014 - ISETBIO', 'VMA mRGCf density Watson 2014 - ISETBIO', ...
                    'HVA Cones Curcio et al 1990 - rgcDisplacementMap', 'VMA Cones Curcio et al 1990 - rgcDisplacementMap', ...        
                    'HVA mRGC Curcio & Allen 1990 - rgcDisplacementMap', 'VMA mRGC Curcio & Allen 1990 - rgcDisplacementMap', ...
                    'HVA V1/V2 cortex', 'VMA V1/V2 cortex'}, 'Location', 'SouthEast');



ylabel('more vertical/upper VF <- Asymmetry (%) -> more horizontal/lower VF')
xlabel('Eccentricity (deg)')
title('Cones vs mRGC RF vs V1 cortex : HVA VMA')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF__vs_V1_CMF'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF__vs_V1_CMF'), '-dpdf', '-fillpage')


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

