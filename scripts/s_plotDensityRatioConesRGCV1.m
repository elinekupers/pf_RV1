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

conesCurcioIsetbio    = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
conesSongIsetbioYoung = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');
conesSongIsetbioOld   = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Old');

for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
end

[fH1, fH2, fH3, fH4] = visualizeConeDensityVsEccenIsetbio(conesCurcioIsetbio, conesSongIsetbioYoung, ...
                       conesSongIsetbioOld, meridianIdx, eccDeg, saveFigures);

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
    load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplacementMaps.mat'), 'allMaps')

else
    % Get data (by computation, takes about 20 minutes)
    [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
        getConeAndRGCDensityDisplacementMap(0.01, 5, 40, ...
            10, true, true);
end

% Plot density along cardinal meridians vs eccen
for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(sampleResPolAng==cardinalMeridianAngles(ii));
end

% pixelsPerDegVisual = 10;
% imRdim = (max(regularSupportPosDegVisual) * pixelsPerDegVisual * 2) -1;
% x = linspace(-1.* (max(regularSupportPosDegVisual)-.1),(max(regularSupportPosDegVisual)-.1), imRdim)
% xy = meshgrid(x);
% 
% figure; surf(x,x,log10(allMaps{3}.data))

[fH5, fH6, fH7, fH8] = visualizeConeAndRGCDensityVsEccenDisplacementToolbox(coneDensityByMeridian, ...
                            mRFDensityByMeridian, meridianIdx, regularSupportPosDegVisual, fH1, fH3, fH4, saveFigures);

                        
                        
%% -----------------------------------------------------------------
%  ------------ CONES, RGCs & V1 from Noah's SFN poster --------------
%  -----------------------------------------------------------------

% Get data
[conesSong, coneDataMeridiansIntegral15] = getConesSong;
[rgcWatson, rgcDataMeridiansIntegral15]  = getMRGCRFWatson;

[fH9, fH10, fH11, fH12] = visualizeConeAndRGCDensitySFNPoster(conesSong, ...
                rgcWatson, coneDataMeridiansIntegral15, rgcDataMeridiansIntegral15, fH7, saveFigures);

% ------------ Plot HVA and VMA points vs eccen for V1 CMF ------------

% Eccentricity range
eccDeg = 0:0.05:40; % deg

% Get CMF from HCP:
v1CMF = getV1CMFHCP;

% Compute HVA and VMA for different wedge sizes
asymV1CMF = struct();

fn = fieldnames(v1CMF.individualSubjects);

for ii = 1:numel(fn)
    
    theseData = v1CMF.individualSubjects.(fn{ii});
    
    horz = mean([theseData(1:181,1), theseData(182:end,1)],2);
    vert = mean([theseData(1:181,2), theseData(182:end,2)],2);
    upr  = mean([theseData(1:181,3), theseData(182:end,3)],2);
    lowr = mean([theseData(1:181,4), theseData(182:end,4)],2);
    
    asymV1CMF.(['hvaAll' fn{ii}]) = 100.* ((horz - vert) ./ mean([horz,vert],2));
    asymV1CMF.(['vmaAll' fn{ii}]) = 100.* ((lowr - upr) ./ mean([upr,lowr],2));
end

% Get CMF from & Hoyt (1991)
CMF_HH91 = HortonHoytCMF(eccDeg);
% Get CMF from Rovamo & Virsu (1979)
CMF_RV79 = RovamuVirsuCMF(eccDeg);
meridCMF_RV79 = [CMF_RV79.nasalR; CMF_RV79.superiorR; CMF_RV79.temporalR; CMF_RV79.inferiorR];


[fH13, fH14, fH15, fH16] = visualizeV1CMFHCP(asymV1CMF, meridCMF_RV79, ...
                CMF_HH91, eccDeg, saveFigures, figureDir);

