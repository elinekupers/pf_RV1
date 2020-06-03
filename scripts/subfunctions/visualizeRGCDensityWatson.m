function [fH9, fH10] = visualizeRGCDensityWatson(rgcWatson, eccDeg, fH8, figureDir, saveFigures)


% ------ Visualize density and HVA vs VMA ------ 
titleStr = 'mRGC RF density Watson 2014 - ISETBIO';
fH9 = plotMeridiansVsEccen(rgcWatson, [], eccDeg, false, titleStr, [], figureDir, saveFigures);


titleStr = 'HVA VMA mRGCf density Watson 2014 - ISETBIO';
fH10 = plotHVAandVMA(rgcWatson, [], eccDeg, false, titleStr, figureDir, saveFigures);

figure(fH10);  % Watson data
h = get(fH8, 'Children'); 

axes(fH10.Children(2))
plot(h(2).Children(2).XData,h(2).Children(2).YData, ...
        'Color', 'b', ...
        'LineWidth', h(2).Children(2).LineWidth, ...
        'LineStyle', h(2).Children(2).LineStyle); hold on
l = findobj(gca, 'Type', 'Line');
legend(l([3,1]), {'HVA mRGCf density Watson 2014 - ISETBIO', ... 
                  'HVA mRGCf Barnett & Aguirre 2018 - rgcDisplacementMap'}, 'Location', 'SouthEast');
legend boxoff;
title('mRGC HVA')
set(gca, 'YLim', [-20, 100]);


axes(fH10.Children(3))
plot(h(1).Children(2).XData,h(1).Children(2).YData, ...
        'Color', 'b', ...
        'LineWidth', h(1).Children(2).LineWidth, ...
        'LineStyle', h(1).Children(2).LineStyle); hold on

    l = findobj(gca, 'Type', 'Line');
    legend(l([3,1]), {'VMA mRGCf density Watson 2014 - ISETBIO', ... 
                      'VMA mRGCf Barnett & Aguirre 2018 - rgcDisplacementMap'}, 'Location', 'SouthEast');
    legend boxoff;
    title('mRGC VMA')
    set(gca, 'YLim', [-10, 80]);

savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_rgcDisplacementMap_Watson'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_rgcDisplacementMap_Watson'), '-dpdf', '-fillpage')
