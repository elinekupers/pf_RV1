%% s_plotDensityConesRGCV1_Benson_NN

% Script to plot cone and RGC density for meridians (either just on the
% meridia or as the integral of a wedge +/- 15 deg polar angle wedge, from
% 1-6 deg eccentricity.


cd(pfRV1rootPath)
saveData               = false;
saveFigures            = true;
loadDataFromServer     = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'NatureNeuro');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end

% Figure options
colors = {'r', 'b', 'g', 'k'};
yl     = [1e2, 3e4]; % y axis limit for cone density in counts/deg2

% Unit converters
deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Meridia angles
cardinalMeridianAngles = [0 90 180 270]; % deg: 0=nasal, 90=superior, 180=temporal and 270=inferior retina of left eye

% Get eccentricity ranges  (1-6 degrees)
eccenBoundary = [0,8];

% Eccentricity range
dt = 0.1;
eccDeg = eccenBoundary(1):dt:eccenBoundary(2); % deg

% Polar angle range
angDeg = 0:0.1:360;  % deg
[~,meridians] = intersect(angDeg,cardinalMeridianAngles);

%% -----------------------------------------------------------------
%  --------- CONES from Curcio et al (1991) using ISETBIO -----------
%  -----------------------------------------------------------------

% Load data from ISETBIO (takes a few minutes)
dataSet = 'Curcio1990';
% coneDensity = getConeDensityIsetbio(angDeg, eccDeg, dataSet);

% or load mat-file
% load(fullfile(pfRV1rootPath, 'external', 'data', 'conesSongISETBIO.mat'),'conesSongIsetbioYoung')
load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'coneDensity')


% Limit data to 1-6 deg eccen
eccenLimits = ((eccenBoundary(1)/dt)+1):1:((eccenBoundary(2)/dt)+1);
meridianData.conesCurcioIsetbio = coneDensity(meridians,eccenLimits)';


% ------ Visualize HVA and VMA ------
titleStr = sprintf('Cone density %s - ISETBIO left eye', dataSet);
visualFieldFlag = false;
fH1 = plotHVAandVMA(meridianData.conesCurcioIsetbio', [], eccDeg, visualFieldFlag, titleStr, figureDir, saveFigures);

%% -----------------------------------------------------------------
%  -------------- mRGC from Watson 2015 using ISETBIO --------------
%  -----------------------------------------------------------------

% Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag
% to true to generate several figures of the Watson 2014 paper
WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);

% Compute total RGC RF density along the superior meridian for a number of eccentricities
meridianLabel = {'nasal meridian','superior meridian','temporal meridian','inferior meridian'};
for ii = 1:length(meridianLabel)
    Watson_mRGCRFDensityPerDeg2(ii,:) = WatsonRGCCalc.midgetRGCRFDensity(eccDeg, meridianLabel{ii}, 'RFs per deg2');
    Watson_totalRGCRFDensityPerDeg2(ii,:) = WatsonRGCCalc.totalRGCRFDensity(eccDeg, meridianLabel{ii}, 'RFs per deg2');
end

% ------ Visualize density and HVA vs VMA ------

titleStr = 'mRGCf density Watson 2014 - ISETBIO';
visualFieldFlag = false;
fH2 = plotHVAandVMA(Watson_mRGCRFDensityPerDeg2, [], eccDeg, visualFieldFlag, titleStr, figureDir, saveFigures);

labels = {{'HVA Cones Curcio et al (1991)', '','HVA mRGC Watson (2014)'}, ...
    {'VMA Cones Curcio et al (1991)','','VMA mRGC Watson (2014)'}};

figure(fH1); h = get(fH2, 'Children');
axes(fH1.Children(2))
plot(h(2).Children(2).XData,h(2).Children(2).YData, ...
            'Color', [0.7 0.7 0.7], ...
            'LineWidth', h(2).Children(2).LineWidth, ...
            'LineStyle', h(2).Children(2).LineStyle); hold on
    legend(labels{1}, 'Location', 'EastOutside'); legend boxoff

axes(fH1.Children(3))
plot(h(1).Children(2).XData,h(1).Children(2).YData, ...
            'Color', [0.7 0.7 0.7], ...
            'LineWidth', h(1).Children(2).LineWidth, ...
            'LineStyle', h(1).Children(2).LineStyle); hold on
    legend(labels{2}, 'Location', 'EastOutside'); legend boxoff 


%% -----------------------------------------------------------------
%  -------------- V1/V2 CMF from HCP Retinotopy dataset ------------
%  -----------------------------------------------------------------

% Get CMF from HCP:
v1CMF = getV1CMFHCP;

% Compute HVA and VMA for different wedge sizes
asymV1CMF = struct();

fn = fieldnames(v1CMF.individualSubjects);

asymPrct = @(x1, x2) 100.*((x1-x2)./nanmean([x1,x2],2));
allData = [];
for ii = 1:numel(fn)
    
    theseData = v1CMF.individualSubjects.(fn{ii});
    allData = cat(3,allData,theseData);
    
    horz = nanmean([theseData(1:181,1), theseData(182:end,1)],2);
    vert = nanmean([theseData(1:181,2), theseData(182:end,2)],2);
    upr  = nanmean([theseData(1:181,3), theseData(182:end,3)],2);
    lowr = nanmean([theseData(1:181,4), theseData(182:end,4)],2);
    
    asymV1CMF.(['hvaAll' fn{ii}]) = 100.* ((horz - vert) ./ nanmean([horz,vert],2));
    asymV1CMF.(['vmaAll' fn{ii}]) = 100.* ((lowr - upr) ./ nanmean([upr,lowr],2));
end

% ------------ Plot HVA and VMA vs eccen for HCP mean subject -----------
eccen = mean([0, 1; ... (0.5 degree)
    1, 2; ... (1.5 degree)
    2, 3; ... (2.5 degree)
    3, 4; ... (3.5 degree)
    4, 5; ... (4.5 degree)
    5, 6; ... (5.5 degree)
    6, 7; ... (6.5 degree)
    7, 8],2); % (7.5 degree)

allEccen = NaN(1,length(eccen)*2);
allEccen(1:2:end) = eccen;
allEccen(2:2:end) = eccen;

fn = fieldnames(asymV1CMF);

lw = 1;
nboot = 1000;

selectGroupHVA = 3:2:length(fn)-4;
selectGroupVMA = 4:2:length(fn)-4;

bootDataHVA = NaN(length(selectGroupHVA), nboot);
bootDataVMA = bootDataHVA;
for f = 1:length(selectGroupHVA)
    bootDataHVA(f,:) = bootstrp(nboot, @(x) nanmean(x), asymV1CMF.(fn{selectGroupHVA(f)}));
end

for f = 1:length(selectGroupVMA)
    bootDataVMA(f,:) = bootstrp(nboot, @(x) nanmean(x), asymV1CMF.(fn{selectGroupVMA(f)}));
end

mdHVA  = nanmedian(bootDataHVA,2);
stdHVA = std(bootDataHVA, [],2,'omitnan');
[pHVA,S1] = polyfit(allEccen(selectGroupHVA),mdHVA',2);
[fitHVA, delta1] = polyval(pHVA,1:6, S1);
Rsquared_HVA = 1 - (S1.normr/norm(mdHVA - mean(mdHVA)))^2;

mdVMA  = nanmedian(bootDataVMA,2);
stdVMA = std(bootDataVMA, [],2,'omitnan');
[pVMA, S2] = polyfit(allEccen(selectGroupVMA),mdVMA',2);
[fitVMA, delta2] = polyval(pVMA,1:6, S2);
Rsquared_VMA = 1 - (S2.normr/norm(mdVMA - mean(mdVMA)))^2;


% Plot HCP integral data points for VMA
figure(fH1); 
figure; set(gcf, 'Position', [418, 69, 1905, 872], 'Color', 'w')
axes(fH1.Children(4));
for ii = 1:length(selectGroupHVA)
    plot(allEccen(selectGroupHVA(ii)),mdHVA(ii), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 12);
    errorbar(allEccen(selectGroupHVA(ii)), mdHVA(ii), stdHVA(ii), 'Color', 'k', 'LineWidth', lw+1);
end
l = findobj(gca, 'Type', 'Line');
legend(l([8,6,2]), {'Cones', 'mRGC', 'V1/V2 cortex'}, 'Location', 'EastOutside'); legend boxoff;
plot(1:6, fitHVA, 'k:', 'LineWidth', lw);
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', [-20,100], 'TickDir', 'out', 'FontSize', 14)
title('HVA')

% Plot HCP integral data points and HVA
axes(fH1.Children(4))
for ii = 1:length(selectGroupVMA)
    plot(allEccen(selectGroupVMA(ii)), mdVMA(ii), 'MarkerFaceColor', 'k', 'Color', 'k', 'Marker', 'o', 'LineWidth', lw, 'MarkerSize', 12);
    errorbar(allEccen(selectGroupVMA(ii)), mdVMA(ii), stdVMA(ii), 'Color', 'k', 'LineWidth', lw+1);
end
l = findobj(gca, 'Type', 'Line');
legend(l([8,6,2]), {'Cones', 'mRGC', 'V1/V2 cortex'}, 'Location', 'EastOutside'); legend boxoff;
plot(1:6, fitVMA, 'k:', 'LineWidth', lw)
ylabel('Asymmetry (%)')
xlabel('Eccentricity (deg)')  
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', [-20,100], 'TickDir', 'out', 'FontSize', 14)
title('VMA')



savefig(fullfile(figureDir, 'Figure3_ConesRGCV1_delta10deg_poly2'))
print(fullfile(figureDir, 'Figure3_ConesRGCV1_delta10deg_poly2'), '-depsc')
