function [fH1, fH2, fH3, fH4] = visualizeConeDensityVsEccenIsetbio(conesSongIsetbioYoung, ...
                            conesSongIsetbioOld, conesSongIsetbio, meridianIdx, eccDeg, saveFigures)

% ------------ Visualize meridan cone density vs eccentricity ------------

yl = [1e2, 3e4];
titleStr = 'Cone density Curcio et al 1990 Group 1- ISETBIO left eye';
meridianData.conesCurcioIsetbioYoung = conesSongIsetbioYoung(meridianIdx,:);
fH1 = plotMeridiansVsEccen(meridianData.conesCurcioIsetbioYoung, eccDeg, titleStr, yl, saveFigures);

titleStr = 'Cone density Curcio et al 1990 Group 2 - ISETBIO left eye';
meridianData.conesSongIsetbioOld = conesSongIsetbioOld(meridianIdx,:);
plotMeridiansVsEccen(meridianData.conesSongIsetbioOld, eccDeg, titleStr, yl, saveFigures);

savefig(fullfile(pfRV1rootPath, 'figures', 'Cone_density_Song_et_al_2011_Group_2-_ISETBIO_left_eye'))
print(fullfile(pfRV1rootPath, 'figures', 'Cone_density_Song_et_al_2011_Group_2-_ISETBIO_left_eye'), '-dpdf', '-fillpage')


% Plot meridian plots
titleStr = 'Cone density Song et al 2011 - ISETBIO left eye';
meridianData.conesSongIsetbio = conesSongIsetbio(meridianIdx,:);
fH2 = plotMeridiansVsEccen(meridianData.conesSongIsetbio, eccDeg, titleStr, yl, saveFigures);


% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Curcio et al 1990 - ISETBIO left eye';
fH3 = plotHVAandVMA(meridianData.conesCurcioIsetbio, eccDeg, titleStr, saveFigures);

titleStr = 'HVA VMA cone density Song et al 2011 - ISETBIO left eye';
fH4 = plotHVAandVMA(meridianData.conesSongIsetbio, eccDeg, titleStr, saveFigures);

titleStr = 'HVA VMA cone density Song et al 2011 Group 2 - ISETBIO left eye';
meridianData.conesSongIsetbioOld = conesSongIsetbioOld(meridianIdx,:);
plotHVAandVMA(meridianData.conesSongIsetbioOld, eccDeg, titleStr, saveFigures);