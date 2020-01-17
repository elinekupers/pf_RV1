function [fH1, fH2, fH3, fH4] = visualizeConeDensityVsEccenIsetbio(conesCurcioIsetbio, conesSongIsetbioYoung, ...
                       conesSongIsetbioOld, meridianIdx, eccDeg, figureDir, saveFigures)
                   
                   
% ------------ Visualize meridan cone density vs eccentricity ------------
yl = [1e2, 3e4];

% Plot meridian cone density plots for Curcio
titleStr = 'Cone density Curcio et al 1990 - ISETBIO left eye';
meridianData.conesCurcioIsetbio = conesCurcioIsetbio(meridianIdx,:);
fH1 = plotMeridiansVsEccen(meridianData.conesCurcioIsetbio, eccDeg, titleStr, yl, figureDir, saveFigures);


% Plot meridian cone density plots for Song 
titleStr = 'Cone density Song et al 2011 Group 1 - ISETBIO left eye';
meridianData.conesSongIsetbioYoung = conesSongIsetbioYoung(meridianIdx,:);
fH2 = plotMeridiansVsEccen(meridianData.conesSongIsetbioYoung, eccDeg, titleStr, yl, figureDir, saveFigures);

titleStr = 'Cone density Song et al 2011 Group 2 - ISETBIO left eye';
meridianData.conesSongIsetbioOld = conesSongIsetbioOld(meridianIdx,:);
plotMeridiansVsEccen(meridianData.conesSongIsetbioOld, eccDeg, titleStr, yl, figureDir, saveFigures);


% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'cone density Curcio et al 1990 - ISETBIO left eye';
fH3 = plotHVAandVMA(meridianData.conesCurcioIsetbio, eccDeg, titleStr, figureDir, saveFigures);

titleStr = 'cone density Song et al 2011 - Group 1 - ISETBIO left eye';
fH4 = plotHVAandVMA(meridianData.conesSongIsetbioYoung, eccDeg, titleStr, figureDir, saveFigures);

titleStr = 'cone density Song et al 2011 Group 2 - ISETBIO left eye';
plotHVAandVMA(meridianData.conesSongIsetbioOld, eccDeg, titleStr, figureDir, saveFigures);