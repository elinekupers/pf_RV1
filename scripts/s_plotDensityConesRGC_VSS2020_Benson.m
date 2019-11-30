%% s_plotDensityConesRGC_VSS2020_Benson

% Script to plot cone and RGC density for meridians (either just on the
% meridia or as the integral of a wedge +/- 15 deg polar angle wedge, from
% 1-6 deg eccentricity.


cd(pfRV1rootPath)
saveData               = false;
saveFigures            = true;
loadDataFromServer     = true;

% Make figure dir if doesnt exist
figureDir = fullfile(pfRV1rootPath, 'figures', 'VSS2020');
if ~exist(figureDir, 'dir'); mkdir(figureDir); end

% Figure options
colors = {'r', 'b', 'g', 'k'};
yl     = [1e2, 3e4]; % y axis limit for cone density in counts/deg2

% Unit converters
deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Meridia angles
cardinalMeridianAngles = [0 90 180 270]; % deg: 0=nasal, 90=superior, 180=temporal and 270=inferior retina of left eye

% Eccentricity range
eccDeg = 0:0.1:40; % deg

% Polar angle range
angDeg = 0:0.1:360;  % deg

% Get eccentricity ranges  (1-6 degrees)
eccenBoundary = [1,6];

[theta, rho] = meshgrid(angDeg, eccDeg);
[X,Y] = pol2cart(theta,rho);

eccen = (eccDeg>=eccenBoundary(1)) & (eccDeg<=eccenBoundary(2)); % deg


% Get polar angle indices (the delta angle for each meridian is 15 deg.)
for ii = 1:length(cardinalMeridianAngles)
    meridians(ii) = find(angDeg==cardinalMeridianAngles(ii));
    
    angleBoundary = cardinalMeridianAngles(ii) + [-15,15];
    if any(angleBoundary<0)
        angleBoundary(angleBoundary<0) = 360 + angleBoundary(angleBoundary<0);
    
        ai1 = (angDeg>=angleBoundary(1) & angDeg<=360);
        ai2 = (angDeg>=0 & angDeg<=angleBoundary(2));
        ai(ii,:) = or(ai1,ai2);
    else   
        ai(ii,:) = (angDeg>=angleBoundary(1) & angDeg<=angleBoundary(2));    
    end 
    
    ma = NaN(size(angDeg));
    me = NaN(size(eccDeg));
    ma(ai(ii,:)) = deg2rad(angDeg(ai(ii,:)));
    ma_single(
    me(eccen) = eccDeg(eccen);

    [mt, mr] = meshgrid(ma,me);
    [mx(:,:,ii), my(:,:,ii)] = pol2cart(mt, mr);
end

clear ai1 ai2 angleBoundary


%% -----------------------------------------------------------------
%  --------- CONES from Song et al. (2011) using ISETBIO -----------
%  -----------------------------------------------------------------

% Load data from ISETBIO (takes a few minutes)
dataSet = 'Song2011Young'; % For Curcio data, set dataSet to 'Curcio1990'
% coneDensity = getConeDensityIsetbio(angDeg, eccDeg, dataSet); 

% or load mat-file
load(fullfile(pfRV1rootPath, 'external', 'data', 'conesSongISETBIO.mat'),'conesSongIsetbioYoung')

% regrid to cartesian space
coneDensity = griddata(theta,rho,conesSongIsetbioYoung', X, Y, 'linear');

% ------ Index data to get polar angle wedges ------ 
mask = zeros(size(coneDensity));

mask.nasal = mask(mx(:,:,ii))


% Get cone density for 15 deg wedge around every meridian
wedge_coneDensity.nasal    = coneDensity.*~isnan(mx(:,:,1)).*~isnan(my(:,:,1));
wedge_coneDensity.superior = coneDensity.*~isnan(mx(:,:,2)).*~isnan(my(:,:,2));
wedge_coneDensity.temporal = coneDensity.*~isnan(mx(:,:,3)).*~isnan(my(:,:,3));
wedge_coneDensity.inferior = coneDensity.*~isnan(mx(:,:,4)).*~isnan(my(:,:,4));

% Or get cone density only the meridia themselves
singleMeridian_coneDensity.nasal    = coneDensity(meridians(1),eccen);
singleMeridian_coneDensity.superior = coneDensity(meridians(2),eccen);
singleMeridian_coneDensity.temporal = coneDensity(meridians(3),eccen);
singleMeridian_coneDensity.inferior = coneDensity(meridians(4),eccen);

% Take the integral of cone density in wedge for every meridian
wedge_coneDensity_integralPA15.nasal    = sum(wedge_coneDensity.nasal(:));
wedge_coneDensity_integralPA15.superior = sum(wedge_coneDensity.superior(:));
wedge_coneDensity_integralPA15.temporal = sum(wedge_coneDensity.temporal(:));
wedge_coneDensity_integralPA15.inferior = sum(wedge_coneDensity.inferior(:));



% Take the integral of cone density along line on every meridia itself
singleMeridian_coneDensity_integralPA15.nasal    = sum(singleMeridian_coneDensity.nasal(:));
singleMeridian_coneDensity_integralPA15.superior = sum(singleMeridian_coneDensity.superior(:));
singleMeridian_coneDensity_integralPA15.temporal = sum(singleMeridian_coneDensity.temporal(:));
singleMeridian_coneDensity_integralPA15.inferior = sum(singleMeridian_coneDensity.inferior(:));

% Concatenate data
coneDataMeridians.wedge = [wedge_coneDensity_integralPA15.nasal, ...
                           wedge_coneDensity_integralPA15.superior, ...
                           wedge_coneDensity_integralPA15.temporal, ...
                           wedge_coneDensity_integralPA15.inferior];
                     
coneDataMeridians.singleMeridian = [singleMeridian_coneDensity_integralPA15.nasal, ...
                           singleMeridian_coneDensity_integralPA15.superior, ...
                           singleMeridian_coneDensity_integralPA15.temporal, ...
                           singleMeridian_coneDensity_integralPA15.inferior];

% Print asymmetries 
fprintf('\nCone density from %s, 1-6 deg eccen, +/- 15 deg polar angle wedge \n', dataSet)
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(coneDataMeridians.wedge))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(coneDataMeridians.wedge))

fprintf('\nCone density from %s, 1-6 deg eccen, only on meridian (no WEDGE) \n', dataSet)
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(coneDataMeridians.singleMeridian))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(coneDataMeridians.singleMeridian))

pointEccen = find(Y(eccen)==mean(eccenBoundary));
hva35_cones = hva([singleMeridian_coneDensity.nasal(pointEccen);singleMeridian_coneDensity.superior(pointEccen);singleMeridian_coneDensity.temporal(pointEccen);singleMeridian_coneDensity.inferior(pointEccen)]);
vma35_cones = vma([singleMeridian_coneDensity.nasal(pointEccen);singleMeridian_coneDensity.superior(pointEccen);singleMeridian_coneDensity.temporal(pointEccen);singleMeridian_coneDensity.inferior(pointEccen)]);

fprintf('\nCone density from %s, 3.5 deg eccen, only on meridian (no WEDGE) \n', dataSet)
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva35_cones)
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma35_cones)

% ------ Visualize HVA and VMA ------ 
% 
% % Plot meridian cone density plots for Song 
% titleStr = sprintf('Cone density %s - ISETBIO left eye',dataSet);
% fH1 = plotMeridiansVsEccen(coneDensity(meridians,:), eccDeg, titleStr, yl, figureDir, saveFigures);
% 
% titleStr = sprintf('HVA VMA cone density %s - ISETBIO left eye', dataSet);
% fH2 = plotHVAandVMA(coneDensity(meridians,:), eccDeg, titleStr, figureDir, saveFigures);
% 
% 
% % Get Song et al. (2011) standard errors (cones/mm2)
% eccentricityMM = [0.18, 0.27, 0.36, 0.45, 0.54, 0.72, 0.90, 1.08, 1.35, 1.62, 1.89, 2.16]; 
% 
% Song2011YoungSEinMM2 = [5.4	2.8	2.9	2.7	2.2	1.8 1.2	1.5	1.5	1.4	0.9	0.6; ... nasal density 10^3 cones/mm^2
%             2.3 1.9 2.0 2.1 1.9 1.4 1.2 0.8 0.8 0.9 0.8 0.5; ... superior density 10^3 cones/mm^2
%             7.5	4.0	2.4	2.0	1.7	1.3 1.5	1.4	0.9	1.0	0.5	0.7; ... temporal density 10^3 cones/mm^2
%             4.5	1.4 3.3	2.8	2.7	2.1 1.6	1.7	1.1	0.9	0.8	0.7]; % inferior density 10^3 cones/mm^2
% 
% % Convert cone density SE to cones/deg2
% eccentricityDegSE = eccentricityMM .* mm2deg;
% Song2011YoungSEinDeg2 = (10^3) .* Song2011YoungSEinMM2./(mm2deg.^2);
% 
% [~, ~, eccentricityDegMn] = intersect(round(eccentricityDegSE,1),round(eccDeg,1));
% 
% % Add errorbars
% figure(fH1); set(gca, 'XLim', [0 7]); hold on;
% for m = 1:4
%     mn = coneDensity(meridians(m),eccentricityDegMn);
%     se = Song2011YoungSEinDeg2(m,:);
%     errorbar(eccentricityDegSE, mn, se, '.','Color', colors{m})
% end
% 
% figure(fH2); set(gca, 'XLim', [0 7]); hold on;
% plot(mean([1,6]),hva(coneDataMeridians.singleMeridian), 'ko', 'MarkerSize',10, 'MarkerFaceColor', 'k');
% plot(mean([1,6]),vma(coneDataMeridians.singleMeridian), 'ko', 'MarkerSize',10);

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

% Interpolate the axes to get 2D density
Watson_mRGCRFDensityPerDeg2(5,:) = Watson_mRGCRFDensityPerDeg2(1,:);
Watson_totalRGCRFDensityPerDeg2(5,:) = Watson_totalRGCRFDensityPerDeg2(1,:);

for aa = 1:numel(eccDeg)
    mRGCfDensity(:,aa) = interp1([cardinalMeridianAngles, 360], Watson_mRGCRFDensityPerDeg2(:,aa), angDeg, 'linear');
    totalRGCfDensity(:,aa) = interp1([cardinalMeridianAngles, 360], Watson_totalRGCRFDensityPerDeg2(:,aa), angDeg, 'linear');
end

% Select mRGC in wedge for every meridian
wedge_mRGCfWatson.nasal    = mRGCfDensity(ai(1,:),eccen);
wedge_mRGCfWatson.superior = mRGCfDensity(ai(2,:),eccen);
wedge_mRGCfWatson.temporal = mRGCfDensity(ai(3,:),eccen);
wedge_mRGCfWatson.inferior = mRGCfDensity(ai(4,:),eccen);

% Select total RGC density in wedge for every meridian
wedge_totalRGCWatson.nasal    = totalRGCfDensity(ai(1,:),eccen);
wedge_totalRGCWatson.superior = totalRGCfDensity(ai(2,:),eccen);
wedge_totalRGCWatson.temporal = totalRGCfDensity(ai(3,:),eccen);
wedge_totalRGCWatson.inferior = totalRGCfDensity(ai(4,:),eccen);

% Calculate mRGC density integral for +/- 15 deg wedge
wedge_mRGCfWatson_integralPA15.nasal    = trapz(trapz(wedge_mRGCfWatson.nasal));
wedge_mRGCfWatson_integralPA15.superior = trapz(trapz(wedge_mRGCfWatson.superior));
wedge_mRGCfWatson_integralPA15.temporal = trapz(trapz(wedge_mRGCfWatson.temporal));
wedge_mRGCfWatson_integralPA15.inferior = trapz(trapz(wedge_mRGCfWatson.inferior));

% Calculate total RGC density integral for +/- 15 deg wedge
wedge_totalRGCWatson_integralPA15.nasal    = trapz(trapz(wedge_totalRGCWatson.nasal));
wedge_totalRGCWatson_integralPA15.superior = trapz(trapz(wedge_totalRGCWatson.superior));
wedge_totalRGCWatson_integralPA15.temporal = trapz(trapz(wedge_totalRGCWatson.temporal));
wedge_totalRGCWatson_integralPA15.inferior = trapz(trapz(wedge_totalRGCWatson.inferior));

% Or take the integral of cone density along just the meridia
singleMeridian_mRGCfWatson.nasal    = mRGCfDensity(meridians(1),eccen);
singleMeridian_mRGCfWatson.superior = mRGCfDensity(meridians(2),eccen);
singleMeridian_mRGCfWatson.temporal = mRGCfDensity(meridians(3),eccen);
singleMeridian_mRGCfWatson.inferior = mRGCfDensity(meridians(4),eccen);

singleMeridian_mRGCfWatson_integralPA15.nasal    = trapz(singleMeridian_mRGCfWatson.nasal);
singleMeridian_mRGCfWatson_integralPA15.superior = trapz(singleMeridian_mRGCfWatson.superior);
singleMeridian_mRGCfWatson_integralPA15.temporal = trapz(singleMeridian_mRGCfWatson.temporal);
singleMeridian_mRGCfWatson_integralPA15.inferior = trapz(singleMeridian_mRGCfWatson.inferior);

% Concatenate data
mrgcDataMeridians.wedge = [wedge_mRGCfWatson_integralPA15.nasal, ...
                           wedge_mRGCfWatson_integralPA15.superior, ...
                           wedge_mRGCfWatson_integralPA15.temporal, ...
                           wedge_mRGCfWatson_integralPA15.inferior];

mrgcDataMeridians.singleMeridian = [singleMeridian_mRGCfWatson_integralPA15.nasal, ...
                          singleMeridian_mRGCfWatson_integralPA15.superior, ...
                          singleMeridian_mRGCfWatson_integralPA15.temporal, ...
                          singleMeridian_mRGCfWatson_integralPA15.inferior];
                       
totalrgcDataMeridians.wedge = [wedge_totalRGCWatson_integralPA15.nasal, ...
                           wedge_totalRGCWatson_integralPA15.superior, ...
                           wedge_totalRGCWatson_integralPA15.temporal, ...
                           wedge_totalRGCWatson_integralPA15.inferior];

fprintf('\nmRGCf density from Watson (2014), 1-6 deg eccen, +/- 15 deg polar angle wedge around meridia\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(mrgcDataMeridians.wedge))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(mrgcDataMeridians.wedge))

fprintf('\nmRGCf density from Watson (2014), 1-6 deg eccen, only on meridian (no WEDGE)\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(mrgcDataMeridians.singleMeridian))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(mrgcDataMeridians.singleMeridian))

fprintf('\nTOTAL RGCf density from Watson (2014), 1-6 deg eccen, +/- 15 deg polar angle wedge around meridia\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(totalrgcDataMeridians.wedge))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(totalrgcDataMeridians.wedge))

pointEccen = find(eccDeg(eccen)==mean(eccenBoundary));
hva35_mrgc = hva([singleMeridian_mRGCfWatson.nasal(pointEccen);singleMeridian_mRGCfWatson.superior(pointEccen);singleMeridian_mRGCfWatson.temporal(pointEccen);singleMeridian_mRGCfWatson.inferior(pointEccen)]);
vma35_mrgc = vma([singleMeridian_mRGCfWatson.nasal(pointEccen);singleMeridian_mRGCfWatson.superior(pointEccen);singleMeridian_mRGCfWatson.temporal(pointEccen);singleMeridian_mRGCfWatson.inferior(pointEccen)]);

fprintf('\nmRGCf density from Watson (2014), 3.5 deg eccen, only on meridian (no WEDGE)\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva35_mrgc)
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma35_mrgc)

ttl = 'Density integral of 15^o wedges, 1-6^o eccen';
hvas = [hva(coneDataMeridians.wedge), hva(mrgcDataMeridians.wedge)];
vmas = [vma(coneDataMeridians.wedge), vma(mrgcDataMeridians.wedge)];

ttl = 'Density integral of meridia only, 1-6^o eccen';
hvas = [hva(coneDataMeridians.singleMeridian), hva(mrgcDataMeridians.singleMeridian)];
vmas = [vma(coneDataMeridians.singleMeridian), vma(mrgcDataMeridians.singleMeridian)];

ttl = 'Density at meridia only, 3.5^o eccen';
hvas = [hva35_cones,hva35_mrgc];
vmas = [vma35_cones,vma35_mrgc];

figure(100); clf; set(gcf, 'Position', [230   273   370   437], 'Color', 'w');
bar([1,1.3], hvas,  0.3, 'EdgeColor','none','facecolor','r', 'FaceAlpha', 0.8, 'LineWidth',3); hold on
bar([1.1,1.4], vmas,  0.3, 'EdgeColor','none','facecolor','b', 'FaceAlpha', 0.8,'LineWidth',3); hold on
title(ttl)
legend({'HVA', 'VMA'}, 'Location', 'SouthEast'); legend boxoff;
set(gca,'xlim',[0.9,1.5],'ylim',[-10,25], 'FontSize', 15);
set(gca,'XTick', [1.05,1.35], 'XTickLabel',{'Cones', 'mRGC'}, 'TickDir', 'out');
ylabel('Asymmetry %'); box off;


% ------ Visualize density and HVA vs VMA ------ 

% % Plot meridian cone density plots for Song 
% titleStr = 'mRGCf density Watson 2014 - ISETBIO';
% fH1 = plotMeridiansVsEccen(Watson_mRGCRFDensityPerDeg2(1:4,:), eccDeg, titleStr, [], figureDir, saveFigures);
% 
% titleStr = 'HVA VMA mRGCf density Watson 2014 - ISETBIO';
% fH2 = plotHVAandVMA(Watson_mRGCRFDensityPerDeg2(1:4,:), eccDeg, titleStr, figureDir, saveFigures);
% figure(fH2); set(gca, 'XLim', [0 40]); hold on;
% plot(mean([1,6]),hva(mrgcDataMeridians.singleMeridian), 'ko', 'MarkerSize',10, 'MarkerFaceColor', 'k')
% plot(mean([1,6]),vma(mrgcDataMeridians.singleMeridian), 'ko', 'MarkerSize',10)


