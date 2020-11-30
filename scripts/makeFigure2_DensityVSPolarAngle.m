function [] = makeFigure2_DensityVSPolarAngle()
% Plot density for cones, mRGC and V1 CMF for each cardinal meridian as a
% function of eccentricity.

cd(pfRV1rootPath)

cardinalMeridianAngles = [0 90 180 270]; % (nasal, superior, temporal and inferior)
colors = {'r', 'b','g', 'k'};
lineStyles = {'-','-','--', '--'};
xScaleDensity = 'log';
yScaleDensity = 'log';

xl = [0.9 40]; % or [0 40];
if strcmp(yScaleDensity, 'log')
    ylR  = 10.^[1 4]; % 3.5 log units for retinal density
    ylV  = 10.^[-0.9 2.1]; % 3.5 log units for cortical surface area
else
    ylR  = [0 25000];   % for retinal density
    ylV  = [0 250];     % for V1 CMF
end

ylR2 = [0 35];  % for transformations
ylV2 = [0 200];  % for transformations

saveData    = false;
saveFigures = true;
loadDataFromServer = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'Figure2');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end


%% ----------------------------------------
%  --------- CONES from ISETBIO -----------
%  ----------------------------------------

% Eccentricity range
dtEcc  = 0.05;      % deg
maxEcc = 40;        % deg
eccDeg = 0:dtEcc:maxEcc; % deg

% Polar angle range
dtAng  = 90;   % deg
maxAng = 360; % deg
angDeg = (0:dtAng:maxAng-1); % deg, (nasal, superior, temporal and inferior)

if loadDataFromServer
    % Get data from server
    if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesCurcioISETBIO.mat'), 'file')
        dataDir = syncDataFromServer();
    end
    load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesCurcioISETBIO.mat'))
    load(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio','conesSongISETBIO.mat'))
    
    for ii = 1:length(cardinalMeridianAngles)
        [~, meridianIdx(ii)] = find(angDeg==cardinalMeridianAngles(ii));
    end
    
    conesCurcioIsetbio = conesCurcioIsetbio(meridianIdx,:);
else
    % Get cone density data from Curcio et al. (1990) and Song et al. (2011)
    conesCurcioIsetbio    = getConeDensityIsetbio(angDeg, eccDeg, 'Curcio1990');
    conesSongIsetbioYoung = getConeDensityIsetbio(angDeg, eccDeg, 'Song2011Young');
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'conesCurcioISETBIO.mat'), 'conesCurcioIsetbio', 'eccDeg', 'angDeg', '-v7.3');
        save(fullfile(pfRV1rootPath, 'external', 'data', 'isetbio', 'conesSongISETBIO.mat'), 'conesSongIsetbioYoung', 'eccDeg', 'angDeg', '-v7.3');
    end
end



%% -----------------------------------------------------------------
%  --------- CONES and RGC from displacement map toolbox -----------
%  -----------------------------------------------------------------

% if loadDataFromServer
%     % Get data from server
%     if ~exist(fullfile(pfRV1rootPath, 'external', 'data', 'coneDensityByMeridian.mat'), 'file')
%         dataDir = syncDataFromServer();
%     end
%     load(fullfile(pfRV1rootPath, 'external', 'data', 'coneDensityByMeridian.mat'),'coneDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
%     load(fullfile(pfRV1rootPath, 'external', 'data', 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')
%     load(fullfile(pfRV1rootPath, 'external', 'data', 'rgcDisplacementMaps.mat'), 'allMaps')
%     
% else
%     % Get data (by computation, takes about 20 minutes)
%     [allMaps, coneDensityByMeridian, rgcDensityByMeridian, regularSupportPosDegVisual] = ...
%         getConeAndRGCDensityDisplacementMap(dtEcc, dtAng, max(eccDeg), ...
%         10, true, true);
% end

% Plot density along cardinal meridians vs eccen
% for ii = 1:length(cardinalMeridianAngles)
%     [~, meridianIdx(ii)] = find(angDeg(1:end-1)==cardinalMeridianAngles(ii));
% end


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

%% -------------------------------------------
%  ----------------- V1 CMF  -----------------
%  -------------------------------------------

% Get CMF from & Hoyt (1991) in mm2/deg2
HH91 = HortonHoytCMF(eccDeg);

% Get CMF from Rovamo & Virsu (1979) in mm2/deg2
% RV79 = RovamuVirsuCMF(eccDeg);
% 
% CMF_RV79(1,:) = RV79.nasalVF;
% CMF_RV79(2,:) = RV79.superiorVF;
% CMF_RV79(3,:) = RV79.temporalVF;
% CMF_RV79(4,:) = RV79.inferiorVF;

loadDataFromServer=false;
if loadDataFromServer
    % Load from server. To compute from scratch: see getV1CMFHCP.m
    load(fullfile(pfRV1rootPath, 'external', 'data', 'V1AreaHCP_sum.mat'),'mdHVA', 'stdHVA', 'mdVMA', 'stdVMA', 'allUpr', 'allLowr', 'allHorz', 'allVert')
else
    % Get CMF data from HCP all 181 subjects, separate for lh and rh:
    v1Area = getV1SurfaceAreaHCP(10);
    numSubjects = length(v1Area.individualSubjects.eccen1_2)/2;
    
    % Get all fieldnames
    fn = fieldnames(v1Area.individualSubjects);
    % predefine variables
    allUpr    = [];
    allLowr   = [];
    allHorz   = [];
    allVert   = [];
    
    % Loop over eccentricities and visual field regions
    for ii = 1:numel(fn) 
        theseData = v1Area.individualSubjects.(fn{ii});
        % Sum across right and left hemispheres
        horz = nansum([theseData(1:numSubjects,1), theseData((numSubjects+1):end,1)],2);
        vert = nansum([theseData(1:numSubjects,2), theseData((numSubjects+1):end,2)],2);
        upr  = nansum([theseData(1:numSubjects,3), theseData((numSubjects+1):end,3)],2);
        lowr = nansum([theseData(1:numSubjects,4), theseData((numSubjects+1):end,4)],2);
        
        allUpr  = [allUpr upr];
        allLowr = [allLowr lowr];
        
        allHorz = [allHorz horz];
        allVert = [allVert vert];
    end
    
    % Scale factor for dorsal/vert
    scaleROI = allVert./(allLowr + allUpr);
    downscaledLowr = scaleROI.*allLowr;
    downscaledUpr = scaleROI.*allUpr;
    
    % Get visual field size for different eccentricity bins
    areaWedge = @(r1,r2,th) (pi*(r1.^2)*(th./360)) - (pi*(r2.^2)*(th./360));
    areaDeg2_40 = [areaWedge(2,1,40), ...
        areaWedge(3,2,40), ...
        areaWedge(4,3,40), ...
        areaWedge(5,4,40), ...
        areaWedge(6,5,40)];
    areaDeg2_20 = [areaWedge(2,1,20), ...
        areaWedge(3,2,20), ...
        areaWedge(4,3,20), ...
        areaWedge(5,4,20), ...
        areaWedge(6,5,20)];
    
    % COMPUTE CMF (surface area / visual area) PER ECCEN BIN
    allUprCMF  = downscaledUpr./areaDeg2_20; 
    allLowrCMF = downscaledLowr./areaDeg2_20; 
    allHorzCMF = allHorz./areaDeg2_40;
    allVertCMF = allVert./areaDeg2_40;
    
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
    
    if saveData
        save(fullfile(pfRV1rootPath, 'external', 'data','benson2020','V1AreaHCP_sum.mat'), ...
            'allUpr','allLowr', 'allHorz','allVert', ...
            'allUprCMF', 'allLowrCMF', 'allHorzCMF', 'allVertCMF',...
            'medianUpperVF', 'medianLowerVF', 'medianHorzVF', 'medianVertVF', ...
            'seUpperVF', 'seLowerVF', 'seHorzVF', 'seVertVF', 'eccenHCP');
    end
end


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

% Get CMF from & Hoyt (1991)
CMF_HH91 = HortonHoytCMF(eccDeg);

% Get CMF from Rovamo & Virsu (1979)
CMF_RV79 = RovamuVirsuCMF(eccDeg);

meridCMF_RV79 = [CMF_RV79.nasalR; ...
    CMF_RV79.superiorR; ...
    CMF_RV79.temporalR; ...
    CMF_RV79.inferiorR];

%% Transformations

% Compute cone : RGC
conesCurcioIsetbio_wHorz(1,:) = nanmean(conesCurcioIsetbio([1,3],:),1);
conesCurcioIsetbio_wHorz(2,:) = conesCurcioIsetbio(2,:);
conesCurcioIsetbio_wHorz(3,:) = conesCurcioIsetbio(4,:);

watsonMRGC_wHorz(1,:) = nanmean(mRGCRFDensityPerDeg2([1,3],:),1);
watsonMRGC_wHorz(2,:) = mRGCRFDensityPerDeg2(2,:);
watsonMRGC_wHorz(3,:) = mRGCRFDensityPerDeg2(4,:);

cone2mRGC        = conesCurcioIsetbio ./ mRGCRFDensityPerDeg2;
cone2mRGC_wHorz  = conesCurcioIsetbio_wHorz ./ watsonMRGC_wHorz;


% rotate/reflect watson rgc data to get visual field coordinates
mRGCWatsonMeridiaIsetbioVF(1,:) = mRGCRFDensityPerDeg2(3,:); % nasal retina -> temporal retina (=nasal VF) 
mRGCWatsonMeridiaIsetbioVF(2,:) = mRGCRFDensityPerDeg2(4,:); % superior retina -> inferior retina (=superior VF)
mRGCWatsonMeridiaIsetbioVF(3,:) = mRGCRFDensityPerDeg2(1,:); % temporal retina -> nasal retina (=temporal VF)
mRGCWatsonMeridiaIsetbioVF(4,:) = mRGCRFDensityPerDeg2(2,:); % inferior retina-> superior retina (= inferior VF)

% Concatenate HCP fits
HCPFits = [fitHorzVFHCP; fitUpperVFHCP; fitLowerVFHCP];

% Subsample x-axis
eccIdx_1_6_Isetbio = (find(eccDeg==1):find(eccDeg==6));

% Merge horizontal axis and concatenate
mRGCWatsonVF_wHorz = [ nanmean([mRGCWatsonMeridiaIsetbioVF(1,eccIdx_1_6_Isetbio); mRGCWatsonMeridiaIsetbioVF(3,eccIdx_1_6_Isetbio)]); ...
                     mRGCWatsonMeridiaIsetbioVF(2,eccIdx_1_6_Isetbio); ... upper VF
                     mRGCWatsonMeridiaIsetbioVF(4,eccIdx_1_6_Isetbio) ];  % lower  VF

% Compute RGC : V1 CMF from HCP ratio
RGCWatson2V1HCP   = mRGCWatsonVF_wHorz ./ HCPFits;

mRGC2V1_Avg = nanmean(mRGCWatsonMeridiaIsetbioVF,1) ./ CMF_HH91;

%% ------------ Visualize meridan cone density vs eccentricity ------------

fH = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [ 305 39 1225 766])

% Cones
subplot(131); hold on;
for ii = 1:4
    plot(eccDeg, conesCurcioIsetbio(ii,:), colors{ii}, 'LineWidth', 2, 'LineStyle', lineStyles{ii});
end
ylabel('Density (counts/deg^2)');
title('Cone density (Curcio et al. 1990)')
set(gca, 'XLim', xl, 'YLim', ylR, 'XScale', xScaleDensity, 'YScale', yScaleDensity, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
set(gca, 'XTick', [1, 10:10:max(xl)], 'XTickLabel', sprintfc('%d',[1, 10:10:max(xl)]))


% Add legend
legend({'Nasal Retina', 'Superior Retina', 'Temporal Retina', 'Inferior Retina'}, ...
     'Location', 'Best'); legend boxoff;

% mRGC
subplot(132); hold on;
for ii = 1:4
    plot(eccDeg, mRGCRFDensityPerDeg2(ii,:), colors{ii}, 'LineWidth', 2, 'LineStyle', lineStyles{ii});
end
ylabel('Density (counts/deg^2)');
title('mRGC RF density (Watson 2015)')
set(gca, 'XLim', xl, 'YLim', ylR,'XScale', xScaleDensity, 'YScale', yScaleDensity, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
set(gca, 'XTick', [1, 10:10:max(xl)], 'XTickLabel', sprintfc('%d',[1, 10:10:max(xl)]))

% Add legend
legend({'Nasal Retina', 'Superior Retina', 'Temporal Retina', 'Inferior Retina'}, ...
     'Location', 'Best'); legend boxoff;

% V1 CMF
subplot(133); cla; hold on;

for k = 1:5 % Plot HCP data
    errorbar(eccenHCP(k), medianHorzVF(k), seHorzVF(k), 'LineWidth',2 , 'color', 'k');
    plot(eccenHCP(k), medianHorzVF(k), 'ro', 'MarkerFaceColor', 'g','MarkerSize', 5, 'LineWidth',2);
    
    errorbar(eccenHCP(k), medianLowerVF(k), seLowerVF(k), 'LineWidth',2, 'color', 'k');
    plot(eccenHCP(k), medianLowerVF(k), 'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth',2);
    
    errorbar(eccenHCP(k), medianUpperVF(k), seUpperVF(k), 'LineWidth',2 , 'color', 'k');
    plot(eccenHCP(k), medianUpperVF(k), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth',2);
end

cla;
% Plot HH
plot(eccDeg,CMF_HH91, 'k:', 'LineWidth',2); % Plot Horton & Hoyt

% Plot fits
plot(x, fitHorzVFHCP, 'r-', 'LineWidth', 2);
plot(x, fitLowerVFHCP, 'b-', 'LineWidth', 2);
plot(x, fitUpperVFHCP, 'k-', 'LineWidth', 2);

% plot average of fits
plot(x, mean([fitHorzVFHCP;fitLowerVFHCP;fitUpperVFHCP]), 'color', [0.3, 0.3, 0.3], 'LineStyle', '--', 'LineWidth', 3);


% Add legend
l = findobj(gca, 'Type', 'Line');
legend(l([4,3,2,1]),{'Horton & Hoyt 91', 'HCP Horizontal VF', 'HCP Lower VF', ...
        'HCP Upper VF', 'Average HCP'}, 'Location', 'Best'); legend boxoff;
% Make figure pretty
set(gca,  'XLim', xl, 'YLim', ylV, 'XScale', xScaleDensity, 'YScale', yScaleDensity, 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
set(gca, 'XTick', [1, 10:10:max(xl)], 'XTickLabel', sprintfc('%d',[1, 10:10:max(xl)]))

ylabel('CMF (mm^2/deg^2)'); xlabel('Eccentricity (deg)');
box off;
title('V1 CMF (H&H, 91; Benson et al. 18)');

if saveFigures
    % Save matlab fig and pdf
    figName = sprintf('Fig2A_DensityPolarAngle_XYScale%s%s_maxEccen%d',xScaleDensity, yScaleDensity,xl(2));
    savefig(fH, fullfile(figureDir, figName))
    hgexport(fH, fullfile(figureDir, figName))
    print(fH, fullfile(figureDir, figName) , '-dpng')
end


%% Transformations cones to mRGC to V1
fH = figure(2); clf; set(gcf, 'Color', 'w', 'Position', [983   301   594   504])
colors_Horz = {'r', 'b', 'k'};

subplot(211); cla; hold on;
for ii = 1:3
    plot(eccDeg, cone2mRGC_wHorz(ii,:), 'color', colors_Horz{ii}, 'LineWidth', 2);
end

xlabel('Eccentricity (deg)');
ylabel('Ratio');
title('Ratio nasal cone density : mRGC RF density')
set(gca, 'XLim', xl, 'YLim', ylR2,'XScale', 'linear', 'YScale', 'linear', 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','off','YMinorGrid','off','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
set(gca, 'XTick', [0:10:max(xl)], 'XTickLabel', sprintfc('%d',[0:10:max(xl)]))

% Add legend
legend({'Horizontal Retina', 'Superior Retina', 'Inferior Retina'}, ...
     'Location', 'Best'); legend boxoff;

subplot(212);cla; hold all;
plot(eccDeg, mRGC2V1_Avg, 'k:', 'LineWidth', 1);
for ii = 1:3
    plot(x, RGCWatson2V1HCP(ii,:), 'color', colors_Horz{ii}, 'LineWidth', 2);
end
plot(x, mean(RGCWatson2V1HCP,1), 'color', [0.3 0.3 0.3], 'LineStyle', '--','lineWidth',3);

xlabel('Eccentricity (deg)');
ylabel('Ratio (counts/mm^2)');
title('Ratio nasal mRGC RF density : V1 CMF')
set(gca, 'XLim', [0 10], 'YLim', ylV2, 'XScale', 'linear', 'YScale', 'linear', 'TickDir', 'out', ...
         'XGrid','on', 'YGrid','on','XMinorGrid','off','YMinorGrid','off','GridAlpha',0.25, ...
         'LineWidth',1,'FontSize',15); box off;
set(gca, 'XTick', [0:2:10], 'XTickLabel', sprintfc('%d',[0:2:10]))

% Add legend
legend({'Horton & Hoyt 1991', 'Horizontal VF', 'HCP Lower VF', ...
        'HCP Upper VF', 'Average HCP'}, 'Location', 'Best'); legend boxoff;


if saveFigures
    % Save matlab fig and pdf
    figName = sprintf('Fig2B_TransformsPolarAngle_XYScale%s%s_maxEccen%d',xScaleDensity, yScaleDensity,xl(2));
    savefig(fH, fullfile(figureDir, figName))
    hgexport(fH, fullfile(figureDir, figName))
    print(fH, fullfile(figureDir, figName) , '-dpng')
end

%% Plot meridian cone density plots for all sources
% titleStr = 'Cone density Curcio et al 1990 - ISETBIO left eye';
% fH1 = plotMeridiansVsEccen(conesCurcioIsetbio(meridianIdx,:), eccDeg, titleStr, yl, figureDir, saveFigures);

% titleStr = 'Cone density Curcio et al 1990 - rgcDisplacement map toolbox';
% fH2 = plotMeridiansVsEccen(coneDensityByMeridian(meridianIdx,:), regularSupportPosDegVisual, titleStr, yl, figureDir, saveFigures);
% 
% titleStr = 'Cone density Song et al 2011 Group 1 - ISETBIO left eye';
% fH3 = plotMeridiansVsEccen(conesSongIsetbioYoung(meridianIdx,:), eccDeg, titleStr, yl, figureDir, saveFigures);

%% ------------ Visualize meridan mRGC density vs eccentricity ------------

% titleStr = 'mRGC RF density Watson 2014 - ISETBIO';
% fH4 = plotMeridiansVsEccen(mRGCRFDensityPerDeg2, eccDeg, titleStr, [], figureDir, saveFigures);

% titleStr = 'mRGC RF density - rgcDisplacement map toolbox';
% fH5 = plotMeridiansVsEccen(mRFDensityByMeridian(meridianIdx,:), regularSupportPosDegVisual, titleStr, [], figureDir, saveFigures);

%% ------------ Visualize V1 CMF vs eccentricity ------------
% fH6 = figure(6); clf; set(gcf, 'Position', [520, 291, 1004, 507]); hold all;

% % Plot Horton & Hoyt
% plot(eccDeg,CMF_HH91, 'k:', 'LineWidth',2);
% 
% % Plot Romavo & Virsu
% for ii = 1:4
%     plot(eccDeg,meridCMF_RV79(ii,:), 'color', colors{ii}, 'LineWidth',2);
% end
% 
% % Plot HCP data
% for k = 1:5
%     errorbar(eccenHCP(k), medianHorzVF(k), seHorzVF(k), 'LineWidth',2 , 'color', 'k');
%     plot(eccenHCP(k), medianHorzVF(k), 'ro', 'MarkerFaceColor', 'g','MarkerSize', 5, 'LineWidth',2);
%     
%     errorbar(eccenHCP(k), medianLowerVF(k), seLowerVF(k), 'LineWidth',2, 'color', 'k');
%     plot(eccenHCP(k), medianLowerVF(k), 'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth',2);
%     
%     errorbar(eccenHCP(k), medianUpperVF(k), seUpperVF(k), 'LineWidth',2 , 'color', 'k');
%     plot(eccenHCP(k), medianUpperVF(k), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth',2);
% end
% 
% % Add legend
% legend({'Horton & Hoyt 91', ...
%     'Romavo & Virsu 87 - Left (HVF)', ...
%     'Romavo & Virsu 87 - Upper (UVF)', ...
%     'Romavo & Virsu 87 - Right (HVF)', ...
%     'Romavo & Virsu 87 - Lower (LVF)', ...
%     '','HCP Horizontal VF', '', 'HCP Lower VF', '' 'HCP Upper VF'}, 'Location', 'BestOutside')
% legend boxoff;
% 
% % Plot fits
% plot(x, fitUpperVFHCP, 'k:', 'LineWidth', 1);
% plot(x, fitHorzVFHCP, 'r:', 'LineWidth', 1);
% plot(x, fitLowerVFHCP, 'b:', 'LineWidth', 1);
% 
% % Make figure pretty
% set(gca, 'XScale','linear','YScale','log', 'XLim', [0 40], 'YLim', [0 100], 'TickDir', 'out', 'FontSize', 12)
% set(gca, 'XTick', [0, 5, 10:10:40], 'XTickLabel', {'0.5', '5', '10', '20', '30', '40'})
% ylabel('CMF (surface area deg^2 / mm^2)'); xlabel('Eccentricity (deg)');
% box off;
% title('V1 CMF - Horton & Hoyt vs Romavo & Virsu vs HCP');


