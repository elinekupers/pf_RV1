% s_plotDensityRatioConesRGCV1

cd(pfRV1rootPath)


%% -------- CONES from ISETBIO ----------

% Eccentricity range
eccDeg = linspace(0,30,73); % deg

% Polar angle range
angDeg = [(0:5:360), 0];  % deg
% angDeg = deg2rad([0, 90, 180, 270, 360]);

plotConeDensityIsetbio(angDeg, eccDeg)

%% -------- RGC from ISETBIO (Not yet implemented) ----------
% fov = 2; % deg
% plotRGCDensityIsetbio(fov)


%% ---------- CONES and RGC from displacement map toolbox ---------

sampleResEccen      = 1; % deg
sampleResPolAng     = 5; % deg
maxEccen            = 40; % deg
pixelsPerDegVisual  = 10; % resolution for plotting (more pixels=higher resolution)

[allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
    plotConeAndRGCDensityDisplacementMap(sampleResEccen, sampleResPolAng, maxEccen, pixelsPerDegVisual);



%% ---------- CONES from Noah's SFN poster ---------
conesSong = plotConesSong;


%% ---------- RGC from Noah's SFN poster ---------
rgcWatson = plotmRGCRFWatson;


%% ---------- V1 CMF from Noah's SFN poster ---------
v1CMF = plotV1CMFHCP;

%% Plot HVA and VMA as a function of eccentricity



