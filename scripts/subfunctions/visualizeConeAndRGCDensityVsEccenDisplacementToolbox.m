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
titleStr = 'HVA VMA cone density Curcio et al 1990 - rgcDisplacement map';
fH7 = plotHVAandVMA(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, figureDir, saveFigures);

titleStr = 'HVA VMA mRGC RF density - rgcDisplacement map';
fH8 = plotHVAandVMA(meridianData.rgcDisplMap, regularSupportPosDegVisual, titleStr, figureDir, saveFigures);


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

% ----- Compare HVA vs VMA vs eccen for cones using ISETBIO vs rgcDisplacement Toolbox --------
figure(fH3);
h2 = get(fH7, 'Children');
h21 = get(fH4, 'Children');  % Song HVA/VMA data

for ii = [3,2]
    plot(h2(2).Children(ii).XData,h2(2).Children(ii).YData, ...
        'Color', 'r', ...
        'LineWidth', h2(2).Children(ii).LineWidth, ...
        'LineStyle', h2(2).Children(ii).LineStyle); hold on
end

legend({'HVA Cones Curcio et al (1990) - Isetbio', 'VMA Cones Curcio et al (1990) - Isetbio', '',...
    'HVA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'VMA Cones Curcio et al (1990) - rgcDisplacement map'}, 'Location', 'SouthEast');

title('ISETBIO (black) vs DisplacementMap (red): Cone HVA VMA Curcio et al 1990')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_left_eye_vs_rgcDisplacementMap'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_left_eye_vs_rgcDisplacementMap'), '-dpdf', '-fillpage')

figure(fH3);
for ii = [3,2]
    plot(h21(2).Children(ii).XData,h21(2).Children(ii).YData, ...
        'Color', 'g', ...
        'LineWidth', h21(2).Children(ii).LineWidth, ...
        'LineStyle', h21(2).Children(ii).LineStyle); hold on
    
end

legend({'HVA Cones Curcio et al (1990) - ISETBIO', 'VMA Cones Curcio et al (1990) - ISETBIO', '',...
    'HVA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'VMA Cones Curcio et al (1990) - rgcDisplacement map', ...
    'HVA Cones Song et al (2011) - ISETBIO', ...
    'VMA Cones Song et al (2011) - ISETBIO'}, 'Location', 'SouthEast');

set(gca, 'YLim', [-30, 30]);
title('Cone density comparisons HVA VMA')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_left_eye_vs_rgcDisplacementMap_vs_Song_et_al_2011'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_left_eye_vs_rgcDisplacementMap_vs_Song_et_al_2011'), '-dpdf', '-fillpage')


% ----- Compare HVA Cones vs mRGC RF vs eccen using rgcDisplacement Toolbox --------
figure(fH7);
h2(2).Children(3).Color = [1 0 0];
h2(2).Children(2).Color = [1 0 0];

h3 = get(fH8, 'Children');
for ii = [3,2]
    plot(h3(2).Children(ii).XData,h3(2).Children(ii).YData, ...
        'Color', 'b', ...
        'LineWidth', h3(2).Children(ii).LineWidth, ...
        'LineStyle', h3(2).Children(ii).LineStyle); hold on
end

legend({'HVA Cones Curcio et al (1990) - rgcDisplacement map', 'VMA Cones Curcio et al (1990) - rgcDisplacement map', '',...
    'HVA mRGC RF - rgcDisplacement map', 'VMA mRGC RF - rgcDisplacement map'}, 'Location', 'SouthEast');


title('Cones (red) vs mRGC RF (blue): HVA VMA rgc displacement map')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'), '-dpdf', '-fillpage')
