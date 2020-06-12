% s_ComputeTransformationRatios

% Script to plot the Cone:mRGC ratio and mRGC:V1 ratio as a function of
% eccentricity.

% Dir to save figures
figureDir = fullfile(pfRV1rootPath, 'figures',  'ratiosConesRGCV1');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end

saveFigs = true;

% Polar angles
meridia = [0 90 180 270]; % deg
angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};
colors = {'r', 'b', 'g', 'k'};

%% -------------------------
% ----- LOAD Cone data ----- 
% --------------------------

% Computation for left eye retina
% eccDeg = 0:0.05:60;
% coneDensity = getConeDensityIsetbio(meridia, eccDeg, 'Curcio1990');

% Curcio et al. 1990 (left eye, retina coords)
load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'conesCurcioIsetbio', 'eccDeg','angDeg');

[~,idxConesLeftEye] = intersect(angDeg, meridia);
idxConesIsetbioRightEye = [idxConesLeftEye(3),idxConesLeftEye(2), idxConesLeftEye(1), idxConesLeftEye(4)];
conesCurcioMeridiaIsetbioRetina = conesCurcioIsetbio(idxConesLeftEye,:);

% rgcDisplacementMap
load(fullfile(pfRV1rootPath, 'external', 'data', 'coneDensityByMeridian.mat'),'coneDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
conesCurcioMeridiaDisplMapRetina = coneDensityByMeridian(idxConesLeftEye,:);

%% -------------------------
% ----- LOAD mRGC data ----- 
% --------------------------

% rgcDisplacementMap
load(fullfile(pfRV1rootPath, 'external', 'data', 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
mRGCMeridiaDisplMapRetina = mRFDensityByMeridian(idxConesLeftEye,:);

% Watson
load(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2')
mRGCWatsonMeridiaIsetbioRetina = mRGCRFDensityPerDeg2([3, 2, 1, 4],:);

% To compute Watson data from scratch:
% Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag
% to true to generate several figures of the Watson 2014 paper
% WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);

% Reproduce meridia data
% meridianLabelWatson          = {'nasal meridian','superior meridian','temporal meridian','inferior meridian'};
% for ii = 1:length(meridianLabelWatson)
%     [~, mRGCWatsonMeridiaIsetbioRetina(ii,:)] = ...
%         WatsonRGCCalc.mRGCRFSpacingAndDensityAlongMeridian(eccDeg, meridianLabelWatson{ii}, ...
%         'deg', 'deg^2', ...
%         'adjustForISETBioConeDensity', false, ...
%         'subtype', 'ONOFF');
% end

% Apparently you can also directly get the cone:RGC ratio.
% Compute total RGC RF density along the superior meridian for a number of eccentricities
% [conesToMRGCratio, spatialSupport, xLabelString, yLabelString, ratioLabel, ...
%     meridianConeToMRGratio, eccUnits] = WatsonRGCCalc.compute2DConeToMRGCRFRatio(eccDeg, 'right eye visual field');
% 
% % transform visual coords to left eye retinal coords
% watsonCone2RGC = [];
% watsonCone2RGC(1,:) = meridianConeToMRGratio.temporal;
% watsonCone2RGC(2,:) = meridianConeToMRGratio.superior';
% watsonCone2RGC(3,:) = meridianConeToMRGratio.nasal;
% watsonCone2RGC(4,:) = meridianConeToMRGratio.inferior';


%% -----------------------------------------
% ----- Compute Cones to mRGC RF RATIO ----- 
% ------------------------------------------

% Ratio for the 4 cardinal meridians in Retina Coordinates
cone2RGC_isetbio  = conesCurcioMeridiaIsetbioRetina ./ mRGCWatsonMeridiaIsetbioRetina;
cone2RGC_displmap = conesCurcioMeridiaDisplMapRetina ./ mRGCMeridiaDisplMapRetina;
% cone2RGC_WatsonDirectlyFromIsetbio   = watsonCone2RGC;

% Or merge nasal and temporal into one horizontal meridian
conesCurcioMeridiaIsetbioRetina_wHorz = conesCurcioMeridiaIsetbioRetina([1,2,3],:);
conesCurcioMeridiaIsetbioRetina_wHorz(1,:) = nanmean(conesCurcioMeridiaIsetbioRetina([1,3],:));

conesCurcioMeridiaDisplMapRetina_wHorz = conesCurcioMeridiaDisplMapRetina([1,2,3],:);
conesCurcioMeridiaDisplMapRetina_wHorz(1,:) = nanmean(conesCurcioMeridiaDisplMapRetina([1,3],:));

mRGCWatsonMeridiaIsetbioRetina_wHorz = mRGCWatsonMeridiaIsetbioRetina([1,2,3],:);
mRGCWatsonMeridiaIsetbioRetina_wHorz(1,:) = nanmean(mRGCWatsonMeridiaIsetbioRetina([1,3],:));

mRGCMeridiaDisplMapRetina_wHorz = mRGCMeridiaDisplMapRetina([1,2,3],:);
mRGCMeridiaDisplMapRetina_wHorz(1,:) = nanmean(mRGCMeridiaDisplMapRetina([1,3],:));

% Compute ratio for 2 vertical and 1 horizontal meridian
cone2RGC_isetbio_wHorz = conesCurcioMeridiaIsetbioRetina_wHorz ./ mRGCWatsonMeridiaIsetbioRetina_wHorz;
cone2RGC_displmap_wHorz = conesCurcioMeridiaDisplMapRetina_wHorz ./ mRGCMeridiaDisplMapRetina_wHorz;


%% --------------------------
% ----- LOAD V1-V2 data ----- 
% ---------------------------

% Get CMF from & Hoyt (1991) in mm2/deg2
CMF_HH91 = HortonHoytCMF(eccDeg);

% Get CMF from Rovamo & Virsu (1979) in mm2/deg2
CMF_RV79 = RovamuVirsuCMF(eccDeg);

meridiaCMF_RV79(1,:) = CMF_RV79.nasalVF;
meridiaCMF_RV79(2,:) = CMF_RV79.superiorVF; 
meridiaCMF_RV79(3,:) = CMF_RV79.temporalVF;
meridiaCMF_RV79(4,:) = CMF_RV79.inferiorVF; 

% Load from server. To compute from scratch: see getV1CMFHCP.m
load(fullfile(pfRV1rootPath, 'external', 'data', 'V1CMF_HCP.mat'),'mdHVA', 'stdHVA', 'mdVMA', 'stdVMA', 'allUpr', 'allLowr', 'allHorz', 'allVert')

% Note 1:
% ventral = upper visual field, corresponding to inferior retina  (4th row)
% dorsal = lower visual field, corresponding to superior retina (2nd row)

% Note 2: 
% Horizontal and vertical surface area are considered the same size, 
% and dorsal/ventral are considered the same size:
% - Horizontal surface area is only in V1.
%   For +/- 10 deg wedge, contains 80 to 90 deg and 90 to 100 (LH), and
%                                  -110 to -80 deg (RH)
% - Vertical surface area is only in V1. 
%   For +/- 10 deg wedge, contains 0 to 10 deg (LH), -10 to 0 deg (RH) and
%                                  170 to 180 deg (LH), -180 to -170 deg (RH)
% - Dorsal and Ventral surface area is V1-V2 combined, containing twice as
% big visual area as vertical. The ventral and dorsal wedges are also 
% combined across LH and RH, but just the upper of lower.
% Since we don't have V2 data, we make the assumption that V2 is roughly
% the same size as V1. So we use V2 to estimate the size of V1. In this
% case we could either scale the ventral and dorsal down to match vertical,
% and divide by the 10 deg (V1-side) wedge only. Or keep the V1-V2 surface
% area and divide by a 20 degree wedge instead.

% Scale factor for dorsal/vert
scaleROI = nanmedian(allVert)./(nanmedian(allLowr) + nanmedian(allUpr));

downscaledLowr = scaleROI.*allLowr;
downscaledUpr = scaleROI.*allUpr;

% Get visual field size for different eccentricity bins
areaWedge = @(r,th) (pi*(r.^2)*(th./360));
areaDeg2 = [areaWedge(2,20) - areaWedge(1,20), ...
    areaWedge(3,20) - areaWedge(2,20), ...
    areaWedge(4,20) - areaWedge(3,20), ...
    areaWedge(5,20) - areaWedge(4,20), ...
    areaWedge(6,20) - areaWedge(5,20)];

% COMPUTE CMF (surface area / visual area) PER ECCEN BIN
% allUprCMF  = downscaledUpr./areaDeg2; % Note if downscaled, we should use 10 deg visual field area, not 20 deg).
% allLowrCMF = downscaledLowr./areaDeg2; % Note if downscaled, we should use 10 deg visual field area, not 20 deg).
allUprCMF  = allUpr./areaDeg2;
allLowrCMF = allLowr./areaDeg2;
allHorzCMF = allHorz./areaDeg2;
allVertCMF = allVert./areaDeg2;

% Boostrap median CMF
nboot = 1000;
bootDataUpperVisualField = bootstrp(nboot, @(x) nanmedian(x), allUprCMF);
bootDataLowerVisualField = bootstrp(nboot, @(x) nanmedian(x), allLowrCMF);
bootDataHorizontalVisualField = bootstrp(nboot, @(x) nanmedian(x), allHorzCMF);
bootDataVerticalVisualField = bootstrp(nboot, @(x) nanmedian(x), allVertCMF);

% Get the median across bootstraps
medianUpperVF = median(bootDataUpperVisualField,1);
medianLowerVF = median(bootDataLowerVisualField,1);
medianHorzVF = median(bootDataHorizontalVisualField,1);
medianVertVF = median(bootDataVerticalVisualField,1);

seUpperVF = std(bootDataUpperVisualField,[],1);
seLowerVF = std(bootDataLowerVisualField,[],1);
seHorzVF = std(bootDataHorizontalVisualField,[],1);
seVertVF = std(bootDataVerticalVisualField,[],1);

% Get eccentricity (to fit data and prepare x-axes for visualization)
eccenHCP = mean([1, 2; ... (1.5 degree)
    2, 3; ... (2.5 degree)
    3, 4; ... (3.5 degree)
    4, 5; ... (4.5 degree)
    5, 6],2); % (5.5 degree)

% Fit a linear function to data
% Note horizontal VF = nasal+temporal retina (1st/3rd row)
x = (1:0.05:6);
[coeffHorzVF]    = fit(eccenHCP,medianHorzVF','power1');
fitHorzVFHCP     = coeffHorzVF.a .* x.^coeffHorzVF.b;

% Note lower VF = superior retina (2th row)
[coeffLowerVF]    = fit(eccenHCP,medianLowerVF','power1');
fitLowerVFHCP     = coeffLowerVF.a .* x.^coeffLowerVF.b;

% Note upper VF = inferior retina (4th row)
[coeffUpperVF]    = fit(eccenHCP,medianUpperVF','power1');
fitUpperVFHCP     = coeffUpperVF.a .* x.^coeffUpperVF.b;



%% DEBUG FIGURE TO COMPARE CMFs of H&H, R&V and HCP data
figure(99); clf; set(gcf, 'Position', [520, 291, 1004, 507]); hold all;

% Plot Horton & Hoyt
plot(eccDeg,CMF_HH91, 'k:', 'LineWidth',2);

% Plot Romavo & Virsu
for ii = 1:4
    plot(eccDeg,meridiaCMF_RV79(ii,:), 'color', colors{ii}, 'LineWidth',2);
end

% Plot HCP data
for k = 1:5
    errorbar(eccenHCP(k), medianHorzVF(k), seHorzVF(k), 'LineWidth',2 , 'color', 'k');
    plot(eccenHCP(k), medianHorzVF(k), 'ro', 'MarkerFaceColor', 'g','MarkerSize', 8, 'LineWidth',2);
    
    errorbar(eccenHCP(k), medianLowerVF(k), seLowerVF(k), 'LineWidth',2, 'color', 'k');
    plot(eccenHCP(k), medianLowerVF(k), 'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth',2);
    
    errorbar(eccenHCP(k), medianUpperVF(k), seUpperVF(k), 'LineWidth',2 , 'color', 'k');
    plot(eccenHCP(k), medianUpperVF(k), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth',2);
end

% Add legend
legend({'Horton & Hoyt 91', ...
    'Romavo & Virsu 87 - Left (HVF)', ...
    'Romavo & Virsu 87 - Upper (UVF)', ...
    'Romavo & Virsu 87 - Right (HVF)', ...
    'Romavo & Virsu 87 - Lower (LVF)', ...
    '','HCP Horizontal VF', '', 'HCP Lower VF', '' 'HCP Upper VF'}, 'Location', 'BestOutside')
legend boxoff;

% Plot fits
plot(x, fitUpperVFHCP, 'k:', 'LineWidth', 1);
plot(x, fitHorzVFHCP, 'r:', 'LineWidth', 1);
plot(x, fitLowerVFHCP, 'b:', 'LineWidth', 1);

% Make figure pretty
set(gca, 'XScale','linear', 'XLim', [0 10], 'YLim', [0 100], 'TickDir', 'out', 'FontSize', 12)
set(gca, 'XTick', [0, 5, 10:10:60], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40', '50','60'})
ylabel('CMF (surface area deg^2 / mm^2)'); xlabel('Eccentricity (deg)');
box off;
title('V1 CMF - Horton & Hoyt vs Romavo & Virsu vs HCP');

if saveFigs
    hgexport(99, fullfile(figureDir,'V1CMF_HCP_HH_VM_linearAxis'));
    print(fullfile(figureDir, 'V1CMF_HCP_HH_VM_linearAxis'), '-dpng');
    savefig(99, fullfile(figureDir,'V1CMF_HCP_HH_VM_linearAxis'));
end

%% -----------------------------------------
% ----- Compute ratio mRGC RF: V1 CMF ------
% ------------------------------------------

% rotate/reflect watson rgc data to get visual field coordinates
mRGCWatsonMeridiaIsetbioVF(1,:) = mRGCWatsonMeridiaIsetbioRetina(1,:); % nasal retina -> temporal retina (=nasal VF) 
mRGCWatsonMeridiaIsetbioVF(2,:) = mRGCWatsonMeridiaIsetbioRetina(4,:); % superior retina -> inferior retina (=superior VF)
mRGCWatsonMeridiaIsetbioVF(3,:) = mRGCWatsonMeridiaIsetbioRetina(3,:); % temporal retina -> nasal retina (=temporal VF)
mRGCWatsonMeridiaIsetbioVF(4,:) = mRGCWatsonMeridiaIsetbioRetina(2,:); % inferior retina-> superior retina (= inferior VF)

% Subsample and rotate/reflect displ map data
subsampleIdx = (1:5:length(regularSupportPosDegVisual));

mRGCMeridiaDisplMapVF(1,:) = mRGCMeridiaDisplMapRetina(1,subsampleIdx);
mRGCMeridiaDisplMapVF(2,:) = mRGCMeridiaDisplMapRetina(4,subsampleIdx);
mRGCMeridiaDisplMapVF(3,:) = mRGCMeridiaDisplMapRetina(3,subsampleIdx);
mRGCMeridiaDisplMapVF(4,:) = mRGCMeridiaDisplMapRetina(2,subsampleIdx);

% Compute ratios mRGC to V1
RGCWatson2V1RV79 =  mRGCWatsonMeridiaIsetbioVF ./ meridiaCMF_RV79;
RGCDisplMap2V1RV79 = mRGCMeridiaDisplMapVF ./ meridiaCMF_RV79;

% Concatenate HCP fits
HCPFits = [fitHorzVFHCP; fitUpperVFHCP; fitLowerVFHCP];

% Subsample x-axis
eccIdx_1_6_Isetbio = (find(eccDeg==1):find(eccDeg==6));
eccDeg_displMap = regularSupportPosDegVisual(subsampleIdx);
eccIdx_1_6_displMap = (find(eccDeg_displMap==1):find(eccDeg_displMap==6));

% Merge horizontal axis and concatenate
mRGCWatsonVF_wHorz = [ nanmean([mRGCWatsonMeridiaIsetbioVF(1,eccIdx_1_6_Isetbio); mRGCWatsonMeridiaIsetbioVF(3,eccIdx_1_6_Isetbio)]); ...
                     mRGCWatsonMeridiaIsetbioVF(2,eccIdx_1_6_Isetbio); ... upper VF
                     mRGCWatsonMeridiaIsetbioVF(4,eccIdx_1_6_Isetbio) ];  % lower  VF

mRGCDisplMapVF_wHorz = [nanmean([mRGCMeridiaDisplMapVF(1,eccIdx_1_6_displMap); mRGCMeridiaDisplMapVF(3,eccIdx_1_6_displMap)]); ...
                     mRGCMeridiaDisplMapVF(2,eccIdx_1_6_displMap); ...  upper VF
                     mRGCMeridiaDisplMapVF(4,eccIdx_1_6_displMap) ]; % lower  VF

% Compute RGC : V1 CMF from HCP ratio
RGCWatson2V1HCP   = mRGCWatsonVF_wHorz ./ HCPFits;
RGCDisplMap2V1HCP = mRGCDisplMapVF_wHorz ./ HCPFits;



%% DEBUG FIGURE TO COMPARE DENSITY TO RATIO

figure(100); clf; set(gcf, 'Position', [520    39   675   759], 'Color', 'w'); 
for ii = 1:3
    subplot(311)
    plot(eccDeg_displMap(eccIdx_1_6_displMap), mRGCDisplMapVF_wHorz(ii,:), 'color', colors{ii}); hold on;
    
    subplot(312)
    plot(eccDeg(eccIdx_1_6_Isetbio), HCPFits(ii,:), 'color', colors{ii}); hold on;
    
     subplot(313)
     ratio = mRGCDisplMapVF_wHorz(ii,:)./HCPFits(ii,:);
    plot(eccDeg(eccIdx_1_6_Isetbio), ratio, 'color', colors{ii}); hold on;
end

subplot(311)
title('mRGC density')
ylabel('mRGC density (cells/deg^2)')
legend({'Horz', 'Upper', 'Lower'})
grid on; set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([0.5 8]); ylim(10.^[0 4])

subplot(312)
title('V1 CMF')
plot(eccDeg,CMF_HH91, 'k:', 'LineWidth',2);
ylabel('V1 CMF (mm^2/deg^2)')
legend({'Horz', 'Upper', 'Lower', 'H&H ''91'}); 
xlim([0.5 8]); ylim(10.^[0 4])
grid on; set(gca, 'XScale', 'log', 'YScale', 'log')

subplot(313)
title('Ratio mRGC:V1')
ylabel('Ratio mRGC:V1 (cells/mm^2)')
legend({'Horz', 'Upper', 'Lower'})
xlabel('Eccentricity (deg)')
xlim([0.5 8]); ylim(10.^[0 4])
grid on; set(gca, 'XScale', 'log', 'YScale', 'log')

if saveFigs
hgexport(100, fullfile(figureDir, 'mRGC2V1CMFRatio_logAxis'));
print(fullfile(figureDir, 'mRGC2V1CMFRatio_logAxis'), '-dpng');
savefig(100, fullfile(figureDir, 'mRGC2V1CMFRatio_logAxis'));
end



%% Visualize  CONE2RGC ratio's

figure(1); set(gcf, 'Position', [961    39   560   759], 'color', 'w'); clf; hold all;
subplot(211)
for ii = 1:4
    plot(eccDeg, cone2RGC_isetbio(ii,:), 'color', colors{ii}, 'lineWidth', 3); hold on;
end
legend(angNames, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 40], 'YLim', [0 30])
% set(gca, 'XTick', [0, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('ISETBIO toolbox - Curcio Cones : Watson mRGC RF');
xlabel('Eccentricity (deg)');
ylabel('Ratio Cone : mRGC RF');
box off; grid on;

subplot(212)
for ii = 1:4
    plot(regularSupportPosDegVisual, cone2RGC_displmap(ii,:), 'color', colors{ii}, 'lineWidth', 3); hold on;
end
legend(angNames, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 40], 'YLim', [0 30])
% set(gca, 'XTick', [0.5, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('Displacement map toolbox - Cone ratio : mRGC RF')
xlabel('Eccentricity (deg)');
ylabel('Ratio Cone : mRGC RF');
box off; grid on;

if saveFigs
    hgexport(1, fullfile(figureDir,'Cone2RGCRatio_linearAxis'));
    print(fullfile(figureDir, 'Cone2RGCRatio_linearAxis'), '-dpng');
    savefig(1, fullfile(figureDir,'Cone2RGCRatio_linearAxis'));
end

%% Plot cone2rgc with averaged horziontal meridian
angNamesHorz = {'Horizontal Retina (HM)', 'Superior Retina (LVM)', 'Inferior Retina (UVM)'};
colorsHorz = {'r', 'b', 'k'};
figure(3); set(gcf, 'Position', [961    39   560   759], 'color', 'w'); clf; hold all;
subplot(211)
for ii = 1:3
    plot(eccDeg, cone2RGC_isetbio_wHorz(ii,:), 'color', colorsHorz{ii}, 'lineWidth', 3); hold on;
end
legend(angNamesHorz, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 40], 'YLim', [0 30])
% set(gca, 'XTick', [0.5, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('ISETBIO toolbox - Curcio Cones : Watson mRGC RF');
xlabel('Eccentricity (deg)');
ylabel('Ratio Cone : mRGC RF');
box off; grid on; axis square

subplot(212)
for ii = 1:3
    plot(regularSupportPosDegVisual, cone2RGC_displmap_wHorz(ii,:), 'color', colorsHorz{ii}, 'lineWidth', 3); hold on;
end
legend(angNamesHorz, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 40], 'YLim', [0 30])
% set(gca, 'XTick', [0.5, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('Displacement map toolbox - Cone ratio : mRGC RF')
xlabel('Eccentricity (deg)');
ylabel('Ratio Cone : mRGC RF');
box off; grid on; axis square

if saveFigs
hgexport(3, fullfile(figureDir, 'Cone2RGCRatio_linearAxis_wHorz'));
print(fullfile(figureDir,'Cone2RGCRatio_linearAxis_wHorz'), '-dpng');
savefig(3,fullfile(figureDir,'Cone2RGCRatio_linearAxis_wHorz'));
end



%% MRGC:V1
figure(2); set(gcf, 'Position', [191 39 1330 759], 'color', 'w'); clf; hold all;
subplot(221)
for ii = 1:4
    plot(eccDeg, RGCWatson2V1RV79(ii,:), 'color', colors{ii}, 'lineWidth', 3); hold on;
end
legend(angNames, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 10], 'YLim', [0 300])
% set(gca, 'XTick', [0.5, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('Watson mRGC RF : Romavo & Virsu 1979');
xlabel('Eccentricity (deg)');
ylabel('Ratio (RGCs/mm^2)');
box off; grid on; axis square

subplot(222)
for ii = 1:4
    plot(eccDeg, RGCDisplMap2V1RV79(ii,:), 'color', colors{ii}, 'lineWidth', 3); hold on;
end
legend(angNames, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 10], 'YLim', [0 300])
% set(gca, 'XTick', [0.5, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
title('Displacement map mRGC RF : Romavo & Virsu 1979')
xlabel('Eccentricity (deg)');
ylabel('Ratio (RGCs/mm^2)');
box off; grid on;  axis square

subplot(223)
subColors = {'r', 'b', 'k'};
for ii = 1:3
    plot(1:0.05:6, RGCWatson2V1HCP(ii,:), 'color', subColors{ii}, 'lineWidth', 3); hold on;
end
legend({'Horizontal VF', 'Upper VF', 'Lower VF'}, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 8], 'YLim', [0 160])
set(gca, 'XTick', [0:2:8], 'XTickLabel', {'0', '2', '4', '6', '8'})
title('HCP 1/x power fit : Watson mRGC RF')
xlabel('Eccentricity (deg)');
ylabel('Ratio (RGCs/mm^2)');
box off; grid on; axis square

subplot(224)
subColors = {'r', 'b', 'k'};
for ii = 1:3
    plot(1:0.05:6, RGCDisplMap2V1HCP(ii,:), 'color', subColors{ii}, 'lineWidth', 3); hold on;
end
legend({'Horizontal VF', 'Upper VF', 'Lower VF'}, 'Location', 'BestOutside'); legend boxoff;
set(gca, 'XScale', 'linear', 'TickDir', 'out', 'FontSize', 12, 'XLim', [0 8], 'YLim', [0 160])
% set(gca, 'XTick', [0:2:8], 'XTickLabel', {'0', '2', '4', '6', '8'})
title('HCP 1/x power fit : Displ mRGC RF ')
xlabel('Eccentricity (deg)');
ylabel('Ratio (RGCs/mm^2)');
box off; grid on; axis square;

if saveFigs
hgexport(2, fullfile(figureDir, 'RGC2V1Ratio_linearAxis'));
print(fullfile(figureDir, 'RGC2V1Ratio_linearAxis'), '-dpng');
savefig(2,fullfile(figureDir,'RGC2V1Ratio_linearAxis'));
end

%% Other rgcDiscplacementMaps for comparison
% load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplacementMaps.mat'), 'allMaps')
