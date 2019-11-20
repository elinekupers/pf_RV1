% s_plotDensityRatioConesRGCV1

cd(pfRV1rootPath)

cardinalMeridianAngles = [0 90 180 270]; % (nasal, superior, temporal and inferior)
saveData    = false;
saveFigures = true;
loadDataFromServer = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end


%% ----------------------------------------
%  --------- CONES from ISETBIO -----------
%  ----------------------------------------

% Eccentricity range
eccDeg = 0:0.05:40; % deg

% Polar angle range
angDeg = [(0:5:360), 0];  % deg

conesCurcioIsetbio = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
conesSongIsetbioYoung   = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');
conesSongIsetbioOld   = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Old');

for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
end

% ------------ Visualize meridan vs eccentricity ------------
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

%% ------------------------------------------------------------
%  --------- (Not yet implemented) RGC from ISETBIO -----------
%  ------------------------------------------------------------

% fov = 2; % deg
% getRGCDensityIsetbio(fov)




%% -----------------------------------------------------------------
%  --------- CONES and RGC from displacement map toolbox -----------
%  -----------------------------------------------------------------

if loadDataFromServer
    % Get data from server
    if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'coneDensityByMeridian.mat'), 'file')
        dataDir = syncDataFromServer();
    end
    load(fullfile(pfRV1rootPath, 'external', 'data', 'coneDensityByMeridian.mat'),'coneDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
    load(fullfile(pfRV1rootPath, 'external', 'data', 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
else
    % Get data (by computation, takes about 20 minutes)
    [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
        getConeAndRGCDensityDisplacementMap();
end

% Plot density along cardinal meridians vs eccen
for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(sampleResPolAng==cardinalMeridianAngles(ii));
end



% ------------ Visualize meridan vs eccentricity ------------
titleStr = 'Cone density Curcio et al 1990 - rgcDisplacement map toolbox';
yl = [1e2, 3e4];
meridianData.conesCurcioDisplMap = coneDensityByMeridian(meridianIdx,:);
fH5 = plotMeridiansVsEccen(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, yl, saveFigures);

titleStr = 'mRGC RF density - rgcDisplacement map toolbox';
meridianData.rgcDisplMap = mRFDensityByMeridian(meridianIdx,:);
fH6 = plotMeridiansVsEccen(meridianData.rgcWatsonDisplMap, regularSupportPosDegVisual, titleStr, [], saveFigures);



% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Curcio et al 1990 - rgcDisplacement map';
fH7 = plotHVAandVMA(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, saveFigures);

titleStr = 'HVA VMA mRGC RF density - rgcDisplacement map';
fH8 = plotHVAandVMA(meridianData.rgcDisplMap, regularSupportPosDegVisual, titleStr, saveFigures);

% ----- Compare Density vs eccen for cones using ISETBIO vs rgcDisplacement Toolbox --------
h1 = get(fH5, 'Children');

figure(fH1);
for ii = 4:-1:1
    plot(h1(2).Children(ii).XData,h1(2).Children(ii).YData, ...
        'Color', h1(2).Children(ii).Color, ...
        'LineWidth', h1(2).Children(ii).LineWidth, ...
        'LineStyle', ':'); hold on
end

title('ISETBIO (solid) vs DisplacementMap (dashed): Cone density Curcio et al 1990')
save(figfullfile(pfRV1rootPath, 'figures', 'Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye')
print(fullfile(pfRV1rootPath, 'figures', 'Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye'), '-dpdf', '-fillpage')

% ----- Compare HVA vs VMA vs eccen for cones using ISETBIO vs rgcDisplacement Toolbox --------
figure(fH3);
h2 = get(fH7, 'Children');
for ii = [3,2]
    plot(h2(2).Children(ii).XData,h2(2).Children(ii).YData, ...
        'Color', [.8 .8 .8], ...
        'LineWidth', h2(2).Children(ii).LineWidth, ...
        'LineStyle', h2(2).Children(ii).LineStyle); hold on
end

title('ISETBIO (black) vs DisplacementMap (red): Cone HVA VMA Curcio et al 1990')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_density_Curcio_et_al_1990_-_ISETBIO_vs_rgcDisplacementMap_left_eye'), '-dpdf', '-fillpage')


% ----- Compare HVA Cones vs mRGC RF vs eccen using rgcDisplacement Toolbox --------
figure(fH7);
h3 = get(fH8, 'Children');
for ii = [3,2]
    plot(h3(2).Children(ii).XData,h3(2).Children(ii).YData, ...
        'Color', 'r', ...
        'LineWidth', h3(2).Children(ii).LineWidth, ...
        'LineStyle', h3(2).Children(ii).LineStyle); hold on
end

title('Cones (black) vs mRGC RF (red): HVA VMA rgc displacement map')
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_Cone_vs_mRGC_RF_density_rgcDisplacementMap'), '-dpdf', '-fillpage')


%% -----------------------------------------------------------------
%  ------------ CONES, RGCs & V1 from Noah's SFN poster --------------
%  -----------------------------------------------------------------

% Get data
[conesSong, coneDataMeridiansIntegral15] = getConesSong;
[rgcWatson, rgcDataMeridiansIntegral15]  = getMRGCRFWatson;

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

% Add integral points
figure(fH11);
hold on; plot(mean([1,6]), hva(coneDataMeridiansIntegral15), 'ro', 'LineWidth', 4);
hold on; plot(mean([1,6]), vma(coneDataMeridiansIntegral15), 'go', 'LineWidth', 4);
legend({'HVA', 'VMA', '', 'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'})
xlim([0, 40])
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_density_Song_et_al_2011_-_SfN_Poster_wIntegral'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_cone_density_Song_et_al_2011_-_SfN_Poster_wIntegral'), '-dpdf', '-fillpage')


figure(fH12);
hold on; plot(mean([1,6]), hva(rgcDataMeridiansIntegral15), 'bo', 'LineWidth', 4);
hold on; plot(mean([1,6]), vma(rgcDataMeridiansIntegral15), 'yo', 'LineWidth', 4);
legend({'HVA', 'VMA', '', 'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'})
xlim([0, 40]);
savefig(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_Watson_2014_-_SfN_Poster_wIntegral'))
print(fullfile(pfRV1rootPath, 'figures', 'HVA_VMA_mRGC_RF_density_Watson_2014_-_SfN_Poster_wIntegral'), '-dpdf', '-fillpage')


%% ------------ Plot HVA and VMA points vs eccen for V1 CMF ------------

% Eccentricity range
eccDeg = 0:0.05:40; % deg

% Get CMF from HCP:
v1CMFMeridiansIntegral15 = getV1CMFHCP;

% Get CMF from & Hoyt (1991)
CMF_HH91 = HortonHoytCMF(eccDeg);
% Get CMF from Romavo & Virsu (1979)
CMF_RV79 = RovamuVirsuCMF(eccDeg);
meridianDataCMF_RV79 = [CMF_RV79.nasalR; CMF_RV79.superiorR; CMF_RV79.temporalR; CMF_RV79.inferiorR];

titleStr = 'CMF Romavo Virsu 1979 - in retinal coords';
fH13 = plotMeridiansVsEccen(meridianDataCMF_RV79, eccDeg, titleStr, [], saveFigures);

figure(fH13); hold on;
plot(eccDeg, CMF_HH91, 'LineWidth', 3, 'Color', [.7 .7 .7])
legend({'nasal meridian on retina',...
    'superior meridian on retina', ...
    'temporal meridian on retina', ...
    'inferior meridian on retina', ...
    'Hoyton & Hort (1991)'}, 'Location', 'Best');

titleStr = 'CMF Romavo Virsu 1979, Hoyton Hort 1991 - in retinal coords';
title(titleStr);
figName = strrep(titleStr,' ','_');
savefig(fH13, fullfile(figureDir, figName))
print(fullfile(figureDir, figName), '-dpdf', '-fillpage')


titleStr = 'HVA VMA CMF Romavo Virsu 1979 - in retinal coords';
fH14 = plotHVAandVMA(meridianDataCMF_RV79, eccDeg, titleStr, saveFigures);


for ii = 1:size(v1CMFMeridiansIntegral15.individualSubjects.eccen0_35,1)
    hvaAll035(ii) = hva(v1CMFMeridiansIntegral15.individualSubjects.eccen0_35(ii,:));
    vmaAll035(ii) = vma(v1CMFMeridiansIntegral15.individualSubjects.eccen0_35(ii,:));
    
    hvaAll16(ii) = hva(v1CMFMeridiansIntegral15.individualSubjects.eccen1_6(ii,:));
    vmaAll16(ii) = vma(v1CMFMeridiansIntegral15.individualSubjects.eccen1_6(ii,:));
    
    hvaAll357(ii) = hva(v1CMFMeridiansIntegral15.individualSubjects.eccen35_7(ii,:));
    vmaAll357(ii) = vma(v1CMFMeridiansIntegral15.individualSubjects.eccen35_7(ii,:));
end

fH = figure(fH14);
eccen = mean([0,3.5; 1,6; 3.5,7],2);

% Plot HCP integral data points
plot(eccen(1), mean(hvaAll035), 'co', 'LineWidth', 4);
plot(eccen(1), mean(vmaAll035), 'mo', 'LineWidth', 4);

plot(eccen(2), mean(hvaAll16), 'co', 'LineWidth', 4);
plot(eccen(2), mean(vmaAll16), 'mo', 'LineWidth', 4);

plot(eccen(3), mean(hvaAll357), 'co', 'LineWidth', 4);
plot(eccen(3), mean(vmaAll357), 'mo', 'LineWidth', 4);

% Plot 0 line
plot(0:8, zeros(1,9), 'k')

xlim([0, 40])
ylim([-80, 80])

ylabel('more vert/inf retina <- Asymmetry (%) -> more horz/sup retina')
xlabel('Eccentricity (deg)')
titleStr = 'HVA VMA V1 CMF HCP mean subject - SfN Poster';
title(titleStr);

legend({'HVA Romavo & Virsu (1979)', ...
    'VMA Romavo & Virsu (1979)', '', ...
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
    savefig(fH, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
    
end



