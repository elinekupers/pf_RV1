function [fH9, fH10] = visualizeRGCDensityWatson(rgcWatson, eccDeg, fH7, figureDir, saveFigures)


% ------ Visualize density and HVA vs VMA ------ 
titleStr = 'mRGC RF density Watson 2014 - ISETBIO';
fH9 = plotMeridiansVsEccen(rgcWatson, eccDeg, titleStr, [], figureDir, saveFigures);


titleStr = 'HVA VMA mRGCf density Watson 2014 - ISETBIO';
fH10 = plotHVAandVMA(rgcWatson, eccDeg, titleStr, figureDir, saveFigures);

figure(fH10); h = get(fH7, 'Children');
for ii = [5,4,2,1]
    plot(h(2).Children(ii).XData,h(2).Children(ii).YData, ...
        'Color', h(2).Children(ii).Color, ...
        'LineWidth', h(2).Children(ii).LineWidth, ...
        'LineStyle', h(2).Children(ii).LineStyle); hold on
end

obj = findobj(gca,'Type','Line');
legend(obj([7,6,4,3,2,1]), {'HVA mRGCf density Watson 2014 - ISETBIO', 'VMA mRGCf density Watson 2014 - ISETBIO', ...
'HVA Cones Curcio et al 1990 - rgcDisplacementMap', 'VMA Cones Curcio et al 1990 - rgcDisplacementMap', ...        
'HVA mRGC Curcio & Allen 1990 - rgcDisplacementMap', 'VMA mRGC Curcio & Allen 1990 - rgcDisplacementMap'}, 'Location', 'SouthEast');

title('Cones vs mRGC RF (Curcio 1990 vs Watson 2014): HVA VMA')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap_Watson'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap_Watson'), '-dpdf', '-fillpage')
