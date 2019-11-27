function [fH9, fH10, fH11, fH12] = visualizeConeAndRGCDensitySFNPoster(conesSong, ...
                rgcWatson, coneDataMeridiansIntegral15, rgcDataMeridiansIntegral15, fH7, saveFigures)

yl = [1e2, 3e4];

            
% ------------ Visualize meridan vs eccentricity ------------
titleStr = 'Cone density Song et al 2011 - SfN Poster';
meridianData.conesSong = [conesSong.nasal; conesSong.superior; conesSong.temporal; conesSong.inferior];
fH9 = plotMeridiansVsEccen(meridianData.conesSong, conesSong.eccentricity, titleStr, yl, saveFigures);

titleStr = 'mRGC RF density Watson 2014 - SfN Poster';
meridianData.rgcWatson = [rgcWatson.nasal; rgcWatson.superior; rgcWatson.temporal; rgcWatson.inferior];
fH10 = plotMeridiansVsEccen(meridianData.rgcWatson, rgcWatson.eccentricity, titleStr, [], saveFigures);


% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Song et al 2011 - SfN Poster';
fH11 = plotHVAandVMA(meridianData.conesSong, conesSong.eccentricity, titleStr, saveFigures);


titleStr = 'HVA VMA mRGC RF density Watson 2014 - SfN Poster';
fH12 = plotHVAandVMA(meridianData.rgcWatson, rgcWatson.eccentricity, titleStr, saveFigures);

figure(fH12);
h7 = get(fH7, 'Children'); hold on;
for ii = [5,4]
    plot(h7(2).Children(ii).XData,h7(2).Children(ii).YData, ...
        'Color', [0 1 0], ...
        'LineWidth', h7(2).Children(ii).LineWidth, ...
        'LineStyle', h7(2).Children(ii).LineStyle); hold on
end

for ii = [2,1]
    plot(h7(2).Children(ii).XData,h7(2).Children(ii).YData, ...
        'Color', h7(2).Children(ii).Color, ...
        'LineWidth', h7(2).Children(ii).LineWidth, ...
        'LineStyle', h7(2).Children(ii).LineStyle); hold on
end
    
xlim([0, 40])
legend({'HVA mRGC RF density Watson 2014 - SfN Poster', 'VMA mRGC RF density Watson 2014 - SfN Poster', '',...
    'HVA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'VMA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'HVA mRGC RF - rgcDisplacement map', ...
    'VMA mRGC RF - rgcDisplacement map'}, 'Location', 'SouthEast');
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_mRGC_density_RF_density_Watson_2014_vs_rgcDisplacement_map_-_SfN_Poster'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_mRGC_density_RF_density_Watson_2014_vs_rgcDisplacement_map_-_SfN_Poster'), '-dpdf', '-fillpage')



% Add integral points
figure(fH11); xlim([0, 40])
hold on; plot(mean([1,6]), hva(coneDataMeridiansIntegral15), 'ro', 'LineWidth', 4);
hold on; plot(mean([1,6]), vma(coneDataMeridiansIntegral15), 'go', 'LineWidth', 4);
legend({'HVA', 'VMA', '', 'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'})
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_density_Song_et_al_2011_-_SfN_Poster_wIntegral'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_density_Song_et_al_2011_-_SfN_Poster_wIntegral'), '-dpdf', '-fillpage')


figure(fH12);
hold on; plot(mean([1,6]), hva(rgcDataMeridiansIntegral15), 'bo', 'LineWidth', 4);
hold on; plot(mean([1,6]), vma(rgcDataMeridiansIntegral15), 'yo', 'LineWidth', 4);
legend({'HVA mRGC RF density Watson 2014 - SfN Poster', 'VMA mRGC RF density Watson 2014 - SfN Poster', '',...
    'HVA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'VMA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'HVA mRGC RF - rgcDisplacement map', ...
    'VMA mRGC RF - rgcDisplacement map', ...
    'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'}, 'Location', 'SouthEast');
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_Watson_2014_vs_rgcDisplacement_map_-_SfN_Poster_wIntegral'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_Watson_2014_vs_rgcDisplacement_map_-_SfN_Poster_wIntegral'), '-dpdf', '-fillpage')
