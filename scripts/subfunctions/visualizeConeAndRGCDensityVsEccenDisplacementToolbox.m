function [fH5, fH6, fH7, fH8] = visualizeConeAndRGCDensityVsEccenDisplacementToolbox(coneDensityByMeridian, ...
                            mRFDensityByMeridian, meridianIdx, regularSupportPosDegVisual, fH1, fH3, fH4, figureDir, saveFigures)

% ------------ Visualize meridan cone density vs eccentricity ------------
titleStr = 'Cone density Curcio et al 1990 - rgcDisplacement map toolbox';
yl = [1e2, 3e4];
meridianData.conesCurcioDisplMap = coneDensityByMeridian(meridianIdx,:);
fH5 = plotMeridiansVsEccen(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, yl, figureDir, saveFigures);

titleStr = 'mRGC RF density - rgcDisplacement map toolbox';
meridianData.rgcDisplMap = mRFDensityByMeridian(meridianIdx,:);
fH6 = plotMeridiansVsEccen(meridianData.rgcDisplMap, regularSupportPosDegVisual, titleStr, [], figureDir, saveFigures);



% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'cone density Curcio et al 1990 - rgcDisplacement map';
fH7  = plotHVAandVMA(meridianData.conesCurcioDisplMap, [], regularSupportPosDegVisual, false, titleStr, figureDir, saveFigures);

titleStr = 'mRGC RF density - rgcDisplacement map';
fH8  = plotHVAandVMA(meridianData.rgcDisplMap, [], regularSupportPosDegVisual,  false, titleStr, figureDir, saveFigures);


% ----- Compare Density vs eccen for cones using ISETBIO vs Song vs rgcDisplacement Toolbox --------
h1 = get(fH5, 'Children'); % curcio density data from displacement map 

figure(fH1);
for ii = 4:-1:1
    plot(h1(2).Children(ii).XData,h1(2).Children(ii).YData, ...
        'Color', h1(2).Children(ii).Color, ...
        'LineWidth', h1(2).Children(ii).LineWidth, ...
        'LineStyle', ':'); hold on
end

title('ISETBIO (solid) vs DisplacementMap (dashed): Cone density Curcio et al 1990')
save(fullfile(pfRV1rootPath, 'figures', 'Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye'))
print(fullfile(pfRV1rootPath, 'figures', 'Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye'), '-dpdf', '-fillpage')

%% ----- Compare HVA vs VMA vs eccen for cones using ISETBIO vs rgcDisplacement Toolbox --------
figure(fH3);  % ISETBIO Curcio HVA data
h2 = get(fH7, 'Children');   % Displacement map Curcio HVA data
h3 = get(fH4, 'Children');  % Song HVA data

ax = get(fH3, 'Children');

axes(fH3.Children(2))
colors = prism(2);
plot(h2(2).Children(2).XData,h2(2).Children(2).YData, ...
        'Color', colors(1,:), ...
        'LineWidth', h2(2).Children(2).LineWidth, ...
        'LineStyle', h2(2).Children(2).LineStyle); hold on
plot(h3(2).Children(2).XData,h3(2).Children(2).YData, ...x
        'Color', colors(2,:), ...
        'LineWidth', h3(2).Children(2).LineWidth, ...
        'LineStyle', h3(2).Children(2).LineStyle); hold on

legend({'HVA Cones Curcio et al. (1990) - Isetbio', '', ... 
        'HVA Cones Curcio et al. (1990) - rgcDispl map', ...
        'HVA Cones Song et al. (2011)   - ISETBIO'}, 'Location', 'SouthEast');
    legend boxoff;

title('Cone density HVA')
set(gca, 'YLim', [-30, 30]);


axes(fH3.Children(3))
plot(h2(1).Children(2).XData,h2(1).Children(2).YData, ...
        'Color', colors(1,:), ...
        'LineWidth', h2(1).Children(2).LineWidth, ...
        'LineStyle', h2(1).Children(2).LineStyle); hold on

plot(h3(1).Children(2).XData,h3(1).Children(2).YData, ...
        'Color', colors(2,:), ...
        'LineWidth', h3(1).Children(2).LineWidth, ...
        'LineStyle', h3(1).Children(2).LineStyle); hold on

legend({'VMA Cones Curcio et al. (1990) - Isetbio', '',...
        'VMA Cones Curcio et al. (1990) - rgcDisplacement map', ...
        'VMA Cones Song et al. (2011)   - ISETBIO'}, 'Location', 'SouthEast'); legend boxoff;
set(gca, 'YLim', [-30, 30]);

savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_Song_et_al_2011_-_ISETBIO_left_eye_vs_rgcDisplacementMap'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_Song_et_al_2011_-_ISETBIO_left_eye_vs_rgcDisplacementMap'), '-dpdf', '-fillpage')



%%

% % ----- Compare HVA Cones vs mRGC RF vs eccen using rgcDisplacement Toolbox --------
% figure(fH7);
% h2(2).Children(3).Color = [1 0 0];
% h2(2).Children(2).Color = [1 0 0];
% 
% h3 = get(fH8, 'Children');
% for ii = [3,2]
%     plot(h3(2).Children(ii).XData,h3(2).Children(ii).YData, ...
%         'Color', 'b', ...
%         'LineWidth', h3(2).Children(ii).LineWidth, ...
%         'LineStyle', h3(2).Children(ii).LineStyle); hold on
% end
% 
% legend({'HVA Cones Curcio et al (1990) - rgcDisplacement map', 'VMA Cones Curcio et al (1990) - rgcDisplacement map', '',...
%     'HVA mRGC RF - rgcDisplacement map', 'VMA mRGC RF - rgcDisplacement map'}, 'Location', 'SouthEast');
% 
% 
% title('Cones (red) vs mRGC RF (blue): HVA VMA rgc displacement map')
% savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'))
% print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'), '-dpdf', '-fillpage')
