function [] = makeSupplementalFig1_HVAVMAvsEccen()

% Radial asymmetries for cone density, mRGC density and V1-V2 surface area
% computed from different publicly available datasets. 
%
% Function was previously a script and called
% s_plotMeridianDensityRatioConesRGCV1vsEccentricity.m

cd(pfRV1rootPath)

cardinalMeridianAngles = [0 90 180 270]; % (nasal, superior, temporal and inferior)
saveData    = false;
saveFigures = false;
loadDataFromServer = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'SupplFig1');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end


%% ----------------------------------------
%  --------- CONES from ISETBIO -----------
%  ----------------------------------------

% Eccentricity range
dtEcc  = 0.05;      % deg
maxEcc = 40;        % deg
eccDeg = 0:dtEcc:maxEcc; % deg

% Polar angle range
dtAng  = 5;   % deg
maxAng = 360; % deg
angDeg = [(0:dtAng:maxAng), 0];  % deg

if loadDataFromServer
    % Get data from server
%     if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesCurcioISETBIO.mat'), 'file')
%         dataDir = syncDataFromServer();
%     end
    load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesCurcioISETBIO.mat'))
    load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesSongISETBIO.mat'))
    
else
    % Get cone density data from Curcio et al. (1990) and Song et al. (2011)
    conesCurcioIsetbio    = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
    conesSongIsetbioYoung = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesCurcioISETBIO.mat'), 'conesCurcioIsetbio', 'eccDeg', 'angDeg', '-v7.3');
        save(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesSongISETBIO.mat'), 'conesSongIsetbioYoung', 'eccDeg', 'angDeg', '-v7.3');
    end
end

% Find cardinal meridian idx if we have fine sampled data
if numel(angDeg)==4
    meridianIdx = 1:4;
else   
    for ii = 1:length(cardinalMeridianAngles)
        [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
    end
end

[fH1, fH2, fH3, fH4] = visualizeConeDensityVsEccenIsetbio(conesCurcioIsetbio, conesSongIsetbioYoung, ...
                       meridianIdx, eccDeg, figureDir, saveFigures);


%% -----------------------------------------------------------------
%  --------- CONES and RGC from displacement map toolbox -----------
%  -----------------------------------------------------------------

if loadDataFromServer
    % Get data from server
    if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplMap', 'coneDensityByMeridian.mat'), 'file')
        dataDir = syncDataFromServer();
    end
    load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplMap', 'coneDensityByMeridian.mat'),'coneDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
    load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplMap', 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
%     load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplMap', 'rgcDisplacementMaps.mat'), 'allMaps')

else
    % Get data (by computation, takes about 20 minutes)
    [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
        getConeAndRGCDensityDisplacementMap(dtEcc, dtAng, max(eccDeg), ...
            10, true, true);
end

% Find cardinal meridian idx if we have fine sampled data
for ii = 1:length(cardinalMeridianAngles)
        [~, meridianIdx_rgcDisp(ii)] = find(sampleResPolAng==cardinalMeridianAngles(ii));
end

% ------------ Plot HVA and VMA vs eccen for Cones and mRGC ------------

[fH5, fH6, fH7, fH8] = visualizeConeAndRGCDensityVsEccenDisplacementToolbox(coneDensityByMeridian, ...
                            mRFDensityByMeridian, meridianIdx_rgcDisp, regularSupportPosDegVisual, fH1, fH3, fH4, figureDir, saveFigures);

%% -----------------------------------------------------------------
%  -------------------- mRGCs from Watson 2014 ---------------------
%  -----------------------------------------------------------------

if loadDataFromServer
    load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2');
else
    % Get mRGC density data per meridian
    mRGCRFDensityPerDeg2  = getMRGCRFWatson(eccDeg);
    
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2', 'eccDeg','cardinalMeridianAngles');
    end
end

% ------------ Plot HVA and VMA vs eccen for mRGC ------------

[fH9, fH10] = visualizeRGCDensityWatson(mRGCRFDensityPerDeg2, eccDeg, fH8, figureDir, saveFigures);

%% -----------------------------------------------------------------
%  ----------------- V1 CMF from Benson et al 2018 -----------------
%  -----------------------------------------------------------------
wwidths = [10,20];
for w = 2:length(wwidths)

% Get surface area data from HCP all 181 subjects, separate for lh and rh:
wedgeWidth = wwidths(w); %  +/- wedge width in deg from meridian
v1CMF = getV1SurfaceAreaHCP(wedgeWidth);

% Get all fieldnames
fn = fieldnames(v1CMF.individualSubjects);

% Get nr of subjects
numSubjects = size(v1CMF.individualSubjects.eccen1_2,2);

% Asymmetry percent diff calculation
asymPrct = @(x1, x2) 100.*((x1-x2)./nanmean([x1,x2],2));

% predefine variables
asymV1CMF = struct();
allUpr    = [];
allLowr   = [];
allHorz   = [];
allVert   = [];

% predefine variables
asymV1CMF = struct();


% Loop over eccentricities and visual field regions
for ii = 1:numel(fn)
    
    % Get subject data for eccen bin
    theseData = v1CMF.individualSubjects.(fn{ii});
    
    % Average across right and left hemispheres
    horz = nanmean(theseData(1,:,:),3);
    vert = nanmean(theseData(2,:,:),3);
    upr  = nanmean(theseData(3,:,:),3);
    lowr = nanmean(theseData(4,:,:),3);
    
    % Compute asymmetry in percent change from mean, for each subject
    asymV1CMF.(['hvaAll' fn{ii}]) = asymPrct(horz,vert);
    asymV1CMF.(['vmaAll' fn{ii}]) = asymPrct(lowr,upr);
    
    allUpr  = [allUpr; upr];
    allLowr = [allLowr; lowr];
    
    allHorz = [allHorz; horz];
    allVert = [allVert; vert];

end

    % Convert zeros to NaN
    allUpr(allUpr==0)=NaN;
    allLowr(allLowr==0)=NaN;
    allHorz(allHorz==0)=NaN;
    allVert(allVert==0)=NaN;

% Get new fieldnames
fn2 = fieldnames(asymV1CMF);

selectGroupHVA = 1:2:length(fn2);
selectGroupVMA = 2:2:length(fn2);

% preallocate space for bootstraps
nboot = 1000;
bootDataHVA = NaN(length(selectGroupHVA), nboot);
bootDataVMA = bootDataHVA;

% Bootstrap the median asymmetry for HVA and VMA across subjects
for f = 1:length(selectGroupHVA)    
    bootDataHVA(f,:) = bootstrp(nboot, @(x) nanmedian(x), asymV1CMF.(fn2{selectGroupHVA(f)}));
end
for f = 1:length(selectGroupVMA)
    bootDataVMA(f,:) = bootstrp(nboot, @(x) nanmedian(x), asymV1CMF.(fn2{selectGroupVMA(f)}));
end

% Compute sample median / se asymmetry across bootstraps, per eccentricity
V1CMF.mdHVA  = median(bootDataHVA,2);
V1CMF.stdHVA = std(bootDataHVA, [],2,'omitnan');

V1CMF.mdVMA  = median(bootDataVMA,2);
V1CMF.stdVMA = std(bootDataVMA, [],2,'omitnan');

if saveData
   save(fullfile(pfRV1rootPath, 'external', 'data', 'benson2021',sprintf('V1CMF_HCP%d.mat', wedgeWidth)),...
       'V1CMF','allUpr', 'allLowr', 'allHorz', 'allVert')
end

% Fit a 2nd degree polynomial to data
% Get eccentricity (to fit data and prepare x-axes for visualization)
V1CMF.polyfit_eccDeg = mean([1, 2; ... (1.5 degree)
    2, 3; ... (2.5 degree)
    3, 4; ... (3.5 degree)
    4, 5; ... (4.5 degree)
    5, 6],2); % (5.5 degree)
V1CMF.polyval_eccDeg = 1:0.05:6; % have a finer sampled eccentricity axis when evaluating fit

[V1CMF.polyfitHVA,V1CMF.polyfitHVAerror]  = polyfit(V1CMF.polyfit_eccDeg,V1CMF.mdHVA,2);
[V1CMF.lineFitHVA, ~] = polyval(V1CMF.polyfitHVA, V1CMF.polyval_eccDeg, V1CMF.polyfitHVAerror);

[V1CMF.polyfitVMA, V1CMF.polyfitVMAerror] = polyfit(V1CMF.polyfit_eccDeg,V1CMF.mdVMA,2);
[V1CMF.lineFitVMA, ~] = polyval(V1CMF.polyfitVMA, V1CMF.polyval_eccDeg, V1CMF.polyfitVMAerror);

% Compute R2
V1CMF.lineFitHVA_R2 = 1 - (V1CMF.polyfitHVAerror.normr / norm(V1CMF.mdHVA - mean(V1CMF.mdHVA)))^2;
V1CMF.lineFitVMA_R2 = 1 - (V1CMF.polyfitVMAerror.normr / norm(V1CMF.mdVMA - mean(V1CMF.mdVMA)))^2;

% Get CMF from & Hoyt (1991)
CMF_HH91 = HortonHoytCMF(eccDeg);

% Get CMF from Rovamo & Virsu (1979)
CMF_RV79 = RovamuVirsuCMF(eccDeg);

meridCMF_RV79 = [CMF_RV79.nasalR; ...
                CMF_RV79.superiorR; ...
                CMF_RV79.temporalR; ...
                CMF_RV79.inferiorR];


% ------------ Plot HVA and VMA points vs eccen for V1 CMF ------------
[fH13, fH14, fH15] = visualizeV1CMFHCP(V1CMF, meridCMF_RV79, ...
                CMF_HH91, eccDeg, saveFigures, figureDir);

end
