% s_plotDensityRatioConesRGCV1

%% Define params

% Eccentricity range
eccDeg = linspace(0,30,73); % deg

% Polar angle range
angDeg = (0:5:360);  % deg
% angDeg = deg2rad([0, 90, 180, 270, 360]);

%% -------- CONES ISETBIO ----------
% plotConeDensityIsetbio(angDeg, eccDeg)

%% -------- RGC ISETBIO (Not yet implemented) ----------
% fov = 2; % deg
% plotConeDensityIsetbio(fov)


%% ---------- CONES, RGC from displacement map toolbox ---------

sampleResEccen      = 1; % deg
sampleResPolAng     = 5; % deg
maxEccen            = 20; % deg
pixelsPerDegVisual  = 10; % resolution for plotting (more pixels=higher resolution)

allMaps = plotConeAndRGCDensityDisplacementMap(sampleResEccen, sampleResPolAng, maxEccen, pixelsPerDegVisual);