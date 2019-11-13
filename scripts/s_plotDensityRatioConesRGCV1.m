% s_plotDensityRatioConesRGCV1

cd(pfRV1rootPath)

cardinalMeridianAngles = [0 90 180 270]; % (nasal, superior, temporal and inferior)
saveData    = false;
saveFigures = true;
loadDataFromServer = true;


%% ----------------------------------------
%  --------- CONES from ISETBIO -----------
%  ----------------------------------------

% Eccentricity range
eccDeg = 0:1:40; % deg

% Polar angle range
angDeg = [(0:5:360), 0];  % deg

conesCurcioIsetbio = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
conesSongIsetbio   = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');

for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
end

% ------------ Visualize meridan vs eccentricity ------------
yl = [1e2, 3e4]; 
titleStr = 'Cone density Curcio Allen 1990 - ISETBIO left eye';
meridianData.conesCurcioIsetbio = conesCurcioIsetbio(meridianIdx,:);
plotMeridiansVsEccen(meridianData.conesCurcioIsetbio, eccDeg, titleStr, yl, saveFigures);

% Plot meridian plots
titleStr = 'Cone density Song et al 2011 - ISETBIO left eye';
meridianData.conesSongIsetbio = conesSongIsetbio(meridianIdx,:);
plotMeridiansVsEccen(meridianData.conesSongIsetbio, eccDeg, titleStr, yl, saveFigures);


% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Curcio Allen 1990 - ISETBIO left eye';
plotHVAandVMA(meridianData.conesCurcioIsetbio, eccDeg, titleStr, saveFigures);

titleStr = 'HVA VMA cone density Song et al 2011 - ISETBIO left eye';
plotHVAandVMA(meridianData.conesSongIsetbio, eccDeg, titleStr, saveFigures);

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
titleStr = 'Cone density Curcio 1990 - rgcDisplacement map';
yl = [1e2, 3e4]; 
meridianData.conesCurcioDisplMap = coneDensityByMeridian(meridianIdx,:);
plotMeridiansVsEccen(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, yl, saveFigures);

titleStr = 'mRGC RF density Watson 2014 - rgcDisplacement map';
meridianData.rgcWatsonDisplMap = mRFDensityByMeridian(meridianIdx,:);
plotMeridiansVsEccen(meridianData.rgcWatsonDisplMap, regularSupportPosDegVisual, titleStr, [], saveFigures);



% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Curcio Allen 1990 - rgcDisplacement map';
plotHVAandVMA(meridianData.conesCurcioDisplMap, regularSupportPosDegVisual, titleStr, saveFigures);

titleStr = 'HVA VMA mRGC RF density Watson 2014 - rgcDisplacement map';
plotHVAandVMA(meridianData.rgcWatsonDisplMap, regularSupportPosDegVisual, titleStr, saveFigures);

%% -----------------------------------------------------------------
%  ------------ CONES, RGCs & V1 from Noah's SFN poster --------------
%  -----------------------------------------------------------------

% Get data 
[conesSong, coneDataMeridiansIntegral15] = getConesSong;
[rgcWatson, rgcDataMeridiansIntegral15]  = getMRGCRFWatson;
v1CMFMeridiansIntegral15 = getV1CMFHCP;

yl = [1e2, 3e4]; 

% ------------ Visualize meridan vs eccentricity ------------
titleStr = 'Cone density Song et al 2011 - SfN Poster';
meridianData.conesSong = [conesSong.nasal; conesSong.superior; conesSong.temporal; conesSong.inferior];
plotMeridiansVsEccen(meridianData.conesSong, conesSong.eccentricity, titleStr, yl, saveFigures);

titleStr = 'mRGC RF density Watson 2014 - SfN Poster';
meridianData.rgcWatson = [rgcWatson.nasal; rgcWatson.superior; rgcWatson.temporal; rgcWatson.inferior];
plotMeridiansVsEccen(meridianData.rgcWatson, rgcWatson.eccentricity, titleStr, [], saveFigures);


% ------------ Plot HVA and VMA vs eccen for cones and mRGC RF ------------
titleStr = 'HVA VMA cone density Song et al 2011 - SfN Poster';
fH1 = plotHVAandVMA(meridianData.conesSong, conesSong.eccentricity, titleStr, saveFigures);

titleStr = 'HVA VMA mRGC RF density Watson 2014 - SfN Poster';
fH2 = plotHVAandVMA(meridianData.rgcWatson, rgcWatson.eccentricity, titleStr, saveFigures);

% Add integral points
figure(fH1);
hold on; plot(mean([1,6]), hva(coneDataMeridiansIntegral15), 'ro', 'LineWidth', 4); 
hold on; plot(mean([1,6]), vma(coneDataMeridiansIntegral15), 'go', 'LineWidth', 4);
legend({'HVA', 'VMA', '', 'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'})
% xlim([0, 40])

figure(fH2);
hold on; plot(mean([1,6]), hva(rgcDataMeridiansIntegral15), 'bo', 'LineWidth', 4); 
hold on; plot(mean([1,6]), vma(rgcDataMeridiansIntegral15), 'yo', 'LineWidth', 4);
legend({'HVA', 'VMA', '', 'HVA integral +/- 15, 1-6 deg eccen', 'VMA integral +/- 15, 1-6 deg eccen'})
% xlim([0, 40])

%% ------------ Plot HVA and VMA points vs eccen for V1 CMF ------------
fH = figure; clf; set(gcf, 'Color', 'w', 'Position', [809   762   751   576]); hold all;
eccen = mean([0,3.5; 1,6; 3.5,7],2);

plot(eccen(1), hva(v1CMFMeridiansIntegral15.eccen0_35), 'co', 'LineWidth', 4); 
plot(eccen(1), vma(v1CMFMeridiansIntegral15.eccen0_35), 'mo', 'LineWidth', 4);

plot(eccen(2), hva(v1CMFMeridiansIntegral15.eccen1_6), 'co', 'LineWidth', 4); 
plot(eccen(2), vma(v1CMFMeridiansIntegral15.eccen1_6), 'mo', 'LineWidth', 4);

plot(eccen(3), hva(v1CMFMeridiansIntegral15.eccen35_7), 'co', 'LineWidth', 4); 
plot(eccen(3), vma(v1CMFMeridiansIntegral15.eccen35_7), 'mo', 'LineWidth', 4);

plot(0:8, zeros(1,9), 'k')

xlim([0, 7])
ylim([-80, 80])

ylabel('more vert/sup retina <- Asymmetry (%) -> more horz/inf retina')
xlabel('Eccentricity (deg)')
titleStr = 'HVA VMA V1 CMF HCP mean subject - SfN Poster';
title(titleStr);

legend({'HVA integral +/- 15, 0-3.5 deg eccen', ...
        'VMA integral +/- 15, 0-3.5 deg eccen', ...
        'HVA integral +/- 15, 1-6 deg eccen', ...
        'VMA integral +/- 15, 1-6 deg eccen', ...
        'HVA integral +/- 15, 3.5-7 deg eccen', ...
        'VMA integral +/- 15, 3.5-7 deg eccen'}, 'Location', 'SouthEast');
legend boxoff;

set(gca,'FontSize', 14,'TickDir', 'out')

if saveFigures
    % Make figure dir if doesnt exist
    figureDir = fullfile(pfRV1rootPath, 'figures');
    if ~exist(figureDir, 'dir'); mkdir(figureDir); end

    % Save matlab fig and pdf
    figName = strrep(titleStr,' ','_');
    savefig(fH, fullfile(figureDir, figName))
    print(fullfile(figureDir, figName), '-dpdf', '-fillpage')
    
end




