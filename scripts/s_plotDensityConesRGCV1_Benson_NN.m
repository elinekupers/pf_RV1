%% s_plotDensityConesRGCV1_Benson_NN
%
% Script to plot cone density, midget RGC density and V1 cortical surface
% area asymmetries along cardinal meridians, as a function of eccentricity.
% Note: Retinal computations are on the meridia, the V1 cortical surface 
% area is the integral of a wedge +/- 10 deg polar angle wedge, from
% 1-6 deg eccentricity.

%% 0. Set up parameters

% Go to rootpath
cd(pfRV1rootPath)

% Define general parameters
saveData               = true;
saveFigures            = true;
loadDataFromServer     = true; % if false, we will recompute the density and surface area numbers (takes 20min+)

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'BensonAsymmetriesNatureNeuro');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end

% Plotting params
yl = [-20,100]; % y-axis limit (percent)
lw = 1;         % line width
visualFieldFlag = false; % false = retinal numbers are in retinal coordinates, not visual field coordinates

% Number of bootstraps for HCP data errorbars
nboot = 1000; 

% Meridia angles
cardinalMeridianAngles = [0 90 180 270]; % deg: 0=nasal, 90=superior, 180=temporal and 270=inferior retina of left eye
meridianLabel          = {'nasal meridian','superior meridian','temporal meridian','inferior meridian'};


%% -----------------------------------------------------------------
%  --------- CONES from Curcio et al (1991) using ISETBIO ----------
%  -----------------------------------------------------------------

dataSet = 'Curcio1990';

% Get eccentricity range
eccenBoundary = [0,8]; % degrees visual angle
dt            = 0.1; % sample rate, degrees visual angle
eccDeg        = eccenBoundary(1):dt:eccenBoundary(2); % degrees visual angle

% Polar angle range
angDeg           = 0:dt:360;  % degrees visual angle
[~,meridiansIdx] = intersect(angDeg,cardinalMeridianAngles);

% Load data from ISETBIO (takes a few minutes)
if loadDataFromServer
    % load data from mat-file
    load(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'coneDensity')
else % recompute data
    coneDensity = getConeDensityIsetbio(angDeg, eccDeg, dataSet);
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'conesCurcioISETBIO.mat'),'coneDensity','eccDeg', 'angDeg','cardinalMeridianAngles','meridianLabel')
    end
end

% Limit data to 1-6 deg eccen for plotting
eccenLimits = ((eccenBoundary(1)/dt)+1):1:((eccenBoundary(2)/dt)+1);
meridianData.conesCurcioIsetbio = coneDensity(meridiansIdx,eccenLimits)';

% ------ Visualize HVA and VMA ------
titleStr = sprintf('Cone density %s - ISETBIO left eye', dataSet);
fH1 = plotHVAandVMA(meridianData.conesCurcioIsetbio', [], eccDeg, visualFieldFlag, titleStr, figureDir, saveFigures);

%% -----------------------------------------------------------------
%  -------------- mRGC from Watson 2015 using ISETBIO --------------
%  -----------------------------------------------------------------

if loadDataFromServer
    % load data from mat-file
    load(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2')
else % recompute data
    % Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag
    % to true to generate several figures of the Watson 2014 paper
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);

    % Compute total RGC RF density along the superior meridian for a number of eccentricities
    for ii = 1:length(meridianLabel)
        [~,mRGCRFDensityPerDeg2(ii,:)] = WatsonRGCCalc.mRGCRFSpacingAndDensityAlongMeridian(eccDeg, meridianLabel{ii}, 'deg', 'deg^2');
    end
    
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'mRGCWatsonISETBIO.mat'),'mRGCRFDensityPerDeg2', 'eccDeg', 'meridianLabel')
    end
end

% ------ Visualize density and HVA vs VMA ------
titleStr = 'mRGCf density Watson 2014 - ISETBIO';
fH2             = plotHVAandVMA(mRGCRFDensityPerDeg2, [], eccDeg, visualFieldFlag, titleStr, figureDir, saveFigures);

% Legend labels
labels = {{'HVA Cones Curcio et al (1990)', '', 'HVA mRGC Watson (2014)'}, ...
          {'VMA Cones Curcio et al (1990)', '', 'VMA mRGC Watson (2014)'}};

% Add mRGC data to cone density figure
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

close(fH2)

%% -----------------------------------------------------------------
%  -------------- V1-V2 CMF from HCP Retinotopy dataset ------------
%  -----------------------------------------------------------------

% Get CMF data from HCP all 181 subjects, separate for lh and rh:
v1CMF = getV1CMFHCP;
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


% Get new fieldnames
fn = fieldnames(asymV1CMF);

selectGroupHVA = 1:2:length(fn);
selectGroupVMA = 2:2:length(fn);

% preallocate space for bootstraps
bootDataHVA = NaN(length(selectGroupHVA), nboot);
bootDataVMA = bootDataHVA;

% Bootstrap the median asymmetry for HVA and VMA across subjects
for f = 1:length(selectGroupHVA)
    bootDataHVA(f,:) = bootstrp(nboot, @(x) nanmedian(x), asymV1CMF.(fn{selectGroupHVA(f)}));
end
for f = 1:length(selectGroupVMA)
    bootDataVMA(f,:) = bootstrp(nboot, @(x) nanmedian(x), asymV1CMF.(fn{selectGroupVMA(f)}));
end

% Compute sample mean / se asymmetry across bootstraps
mdHVA  = mean(bootDataHVA,2);
stdHVA = std(bootDataHVA, [],2,'omitnan');

mdVMA  = mean(bootDataVMA,2);
stdVMA = std(bootDataVMA, [],2,'omitnan');

% Fit a 2nd degree polynomial to data
[pHVA,error1]    = polyfit(allEccen(selectGroupHVA),mdHVA',2);
[fitHVA, delta1] = polyval(pHVA,1:6, error1);

[pVMA, error2]   = polyfit(allEccen(selectGroupVMA),mdVMA',2);
[fitVMA, delta2] = polyval(pVMA,1:6, error2);

% Compute R2
R2_HVA = 1 - (error1.normr/norm(mdHVA - mean(mdHVA)))^2;
R2_VMA = 1 - (error2.normr/norm(mdVMA - mean(mdVMA)))^2;


%  ------------ Plot HVA and VMA vs eccen for HCP mean subject -----------

% prepare x-axes
eccen = mean([1, 2; ... (1.5 degree)
    2, 3; ... (2.5 degree)
    3, 4; ... (3.5 degree)
    4, 5; ... (4.5 degree)
    5, 6],2); % (5.5 degree)

% Duplicate eccentricity x axis, alternate for HVA and VMA
allEccen = NaN(1,length(eccen)*2);
allEccen(1:2:end) = eccen;
allEccen(2:2:end) = eccen;

figure(fH1); 
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
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', yl, 'TickDir', 'out', 'FontSize', 14)
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
set(gca', 'xlim', [0 max(eccDeg)], 'ylim', yl, 'TickDir', 'out', 'FontSize', 14)
title('VMA')


if saveFigures
    savefig(fullfile(figureDir, 'Figure3_ConesRGCV1_delta10deg_poly2'))
    print(fullfile(figureDir, 'Figure3_ConesRGCV1_delta10deg_poly2'), '-depsc')
    print(fullfile(figureDir, 'Figure3_ConesRGCV1_delta10deg_poly2'), '-dpng')
end