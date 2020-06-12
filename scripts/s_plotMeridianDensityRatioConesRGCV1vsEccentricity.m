% s_plotMeridianDensityRatioConesRGCV1vsEccentricity

cd(pfRV1rootPath)

cardinalMeridianAngles = [0 90 180 270]; % (nasal, superior, temporal and inferior)
saveData    = true;
saveFigures = true;
loadDataFromServer = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'DissertationChapter');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end


%% ----------------------------------------
%  --------- CONES from ISETBIO -----------
%  ----------------------------------------

% Eccentricity range
dtEcc  = 0.01;      % deg
maxEcc = 60;        % deg
eccDeg = 0:dtEcc:maxEcc; % deg

% Polar angle range
dtAng  = 5;   % deg
maxAng = 360; % deg
angDeg = [(0:dtAng:maxAng), 0];  % deg

if loadDataFromServer
    % Get data from server
    if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'), 'file')
        dataDir = syncDataFromServer();
    end
    load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'))
    load(fullfile(pfRV1rootPath, 'external', 'data', 'conesSongISETBIO.mat'))
else
    % Get cone density data from Curcio et al. (1990) and Song et al. (2011)
    conesCurcioIsetbio    = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
    conesSongIsetbioYoung = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'), 'conesCurcioIsetbio', 'eccDeg', 'angDeg', '-v7.3');
        save(fullfile(pfRV1rootPath, 'external', 'data', 'conesSongISETBIO.mat'), 'conesSongIsetbioYoung', 'eccDeg', 'angDeg', '-v7.3');
    end
end


for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
end

[fH1, fH2, fH3, fH4] = visualizeConeDensityVsEccenIsetbio(conesCurcioIsetbio, conesSongIsetbioYoung, ...
                       meridianIdx, eccDeg, figureDir, saveFigures);


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
        getConeAndRGCDensityDisplacementMap(dtEcc, dtAng, max(eccDeg), ...
            10, true, true);
end

% Plot density along cardinal meridians vs eccen
for ii = 1:length(cardinalMeridianAngles)
    [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
end

% ------------ Plot HVA and VMA vs eccen for Cones and mRGC ------------

[fH5, fH6, fH7, fH8] = visualizeConeAndRGCDensityVsEccenDisplacementToolbox(coneDensityByMeridian, ...
                            mRFDensityByMeridian, meridianIdx, regularSupportPosDegVisual, fH1, fH3, fH4, figureDir, saveFigures);
               
%% -----------------------------------------------------------------
%  -------------------- mRGCs from Watson 2014 ---------------------
%  -----------------------------------------------------------------

if loadDataFromServer
    load(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2');
else
    % Get mRGC density data per meridian
    mRGCRFDensityPerDeg2  = getMRGCRFWatson(eccDeg);
    
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2', 'eccDeg','cardinalMeridianAngles');
    end
end

% ------------ Plot HVA and VMA vs eccen for mRGC ------------

[fH9, fH10] = visualizeRGCDensityWatson(mRGCRFDensityPerDeg2, eccDeg, fH8, figureDir, saveFigures);

%% -----------------------------------------------------------------
%  ----------------- V1 CMF from Benson et al 2018 -----------------
%  -----------------------------------------------------------------

% Get CMF data from HCP all 181 subjects, separate for lh and rh:
v1CMF = getV1CMFHCP(10);
numSubjects = 181;

% Get all fieldnames
fn = fieldnames(v1CMF.individualSubjects);

% Asymmetry percent diff calculation
asymPrct = @(x1, x2) 100.*((x1-x2)./nanmean([x1,x2],2));

% predefine variables
asymV1CMF = struct();
allData   = [];

% Loop over eccentricities and visual field regions
for ii = 1:numel(fn)
    
    theseData = v1CMF.individualSubjects.(fn{ii});
    allData = cat(3,allData,theseData);
    
    % Average across right and left hemispheres
    horz = nanmean([theseData(1:numSubjects,1), theseData((numSubjects+1):end,1)],2);
    vert = nanmean([theseData(1:numSubjects,2), theseData((numSubjects+1):end,2)],2);
    upr  = nanmean([theseData(1:numSubjects,3), theseData((numSubjects+1):end,3)],2);
    lowr = nanmean([theseData(1:numSubjects,4), theseData((numSubjects+1):end,4)],2);
    
    % Compute asymmetry in percent change from mean, for each subject
    asymV1CMF.(['hvaAll' fn{ii}]) = asymPrct(horz,vert);
    asymV1CMF.(['vmaAll' fn{ii}]) = asymPrct(lowr,upr);
end

if saveData
    save(fullfile(pfRV1rootPath, 'external', 'data', 'V1CMFHPC.mat'), 'asymV1CMF')
end

nboot = 1000; 


% Get CMF from & Hoyt (1991)
CMF_HH91 = HortonHoytCMF(eccDeg);

% Get CMF from Rovamo & Virsu (1979)
CMF_RV79 = RovamuVirsuCMF(eccDeg);

meridCMF_RV79 = [CMF_RV79.nasalR; ...
                CMF_RV79.superiorR; ...
                CMF_RV79.temporalR; ...
                CMF_RV79.inferiorR];


% ------------ Plot HVA and VMA points vs eccen for V1 CMF ------------
[fH13, fH14, fH15] = visualizeV1CMFHCP(asymV1CMF, meridCMF_RV79, ...
                CMF_HH91, eccDeg, fH10, saveFigures, figureDir);

