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

% Get eccentricity ranges  (1-6 degrees)
eccenBoundary = [1,6];

% Eccentricity range
dt = 0.1;
eccDeg = eccenBoundary(1):dt:eccenBoundary(2); % deg

% Polar angle range
angDeg = 0:0.1:360;  % deg

% Make polar angle grid, convert to cartesian grid
[theta, rho] = meshgrid(angDeg, eccDeg);
[X,Y] = pol2cart(theta,rho);


% Get polar angle indices (the delta angle for each meridian is 15 deg.)
for ii = 1:length(cardinalMeridianAngles)
    meridians(ii) = find(angDeg==cardinalMeridianAngles(ii));
    
    angleBoundary = cardinalMeridianAngles(ii) + [-15,15];
    if any(angleBoundary<0)
        angleBoundary(angleBoundary<0) = 360 + angleBoundary(angleBoundary<0);
    
        ai1 = (angDeg>=angleBoundary(1) & angDeg<=360);
        ai2 = (angDeg>0 & angDeg<=angleBoundary(2));
        ai = or(ai1,ai2);
    else   
        ai = (angDeg>=angleBoundary(1) & angDeg<=angleBoundary(2));    
    end 

    % Get polar angle coords for wedge, and convert to cartesian
    angMask = NaN(size(angDeg));
    angMask(ai) = angDeg(ai);
    [thetaM, rhoM] = meshgrid(angMask, eccDeg);
     
    [wedgeMaskX(:,:,ii), wedgeMaskY(:,:,ii)]  = pol2cart(deg2rad(thetaM),rhoM);
    
    % Get polar angle coords for single meridian, and convert to cartesian
    angMask_single = NaN(size(angDeg));
    angMask_single(meridians(ii)) = angDeg(meridians(ii));
    [thetaM_single, rhoM_single] = meshgrid(angMask_single, eccDeg);
    
    [singleMeridianX(:,:,ii),singleMeridianY(:,:,ii)] = pol2cart(deg2rad(thetaM_single),rhoM_single);
    
end

clear ai1 ai2 angleBoundary

figure; 
for ii = 1:4
    subplot(1,2,1); hold all;
    plot(wedgeMaskX(:,:,ii), wedgeMaskY(:,:,ii), 'Color', colors{ii});
    set(gca, 'XLim', [-6 6], 'YLim', [-6 6], 'TickDir', 'out', 'FontSize',12, ...
        'XTick', [-6:6],  'YTick', [-6:6]);  grid on; axis square;
     
    subplot(1,2,2); hold all;
    plot(singleMeridianX(:,:,ii), singleMeridianY(:,:,ii), colors{ii}); 
    set(gca, 'XLim', [-6 6], 'YLim', [-6 6], 'TickDir', 'out', 'FontSize',12, ...
        'XTick', [-6:6],  'YTick', [-6:6]);  grid on; axis square;
end



%% -----------------------------------------------------------------
%  --------- CONES from Song et al. (2011) using ISETBIO -----------
%  -----------------------------------------------------------------

% Load data from ISETBIO (takes a few minutes)
dataSet = 'Song2011Young'; % For Curcio data, set dataSet to 'Curcio1990'
% coneDensity = getConeDensityIsetbio(angDeg, eccDeg, dataSet); 

% or load mat-file
load(fullfile(pfRV1rootPath, 'external', 'data', 'conesSongISETBIO.mat'),'conesSongIsetbioYoung')

% Limit data to 1-6 deg eccen
eccenLimits = ((eccenBoundary(1)/dt)+1):1:((eccenBoundary(2)/dt)+1);
coneDensity = conesSongIsetbioYoung(:,eccenLimits)';

% regrid to polar space
% [Xq, Yq, Vq] = griddata(X,Y,V, xq, yq)
% [X, Y, coneDensity] = griddata(theta,rho,coneDensity, X, Y, 'linear');

% coneDensity = imrotate(coneDensity', -90);

% ------ Index data to get polar angle wedges ------ 
for ii = 1:numel(meridians)
    
    thisWedgeX = wedgeMaskX(:,:,ii);
    thisWedgeY = wedgeMaskY(:,:,ii);
    
    thisSingleMeridianX = singleMeridianX(:,:,ii);
    thisSingleMeridianY = singleMeridianY(:,:,ii);

    wedge_coneDensity(ii,:) = coneDensity(~isnan(thisWedgeX));
    singleMeridian_coneDensity(ii,:) = coneDensity(~isnan(thisSingleMeridianX));
end

% Take the integral of cone density in wedge for every meridian
wedge_coneDensity_integralPA15.nasal    = sum(wedge_coneDensity(1,:),'omitnan');
wedge_coneDensity_integralPA15.superior = sum(wedge_coneDensity(2,:),'omitnan');
wedge_coneDensity_integralPA15.temporal = sum(wedge_coneDensity(3,:),'omitnan');
wedge_coneDensity_integralPA15.inferior = sum(wedge_coneDensity(4,:),'omitnan');

% Take the integral of cone density along line on every meridia itself
singleMeridian_coneDensity_integralPA15.nasal    = sum(singleMeridian_coneDensity(1,:),'omitnan');
singleMeridian_coneDensity_integralPA15.superior = sum(singleMeridian_coneDensity(2,:),'omitnan');
singleMeridian_coneDensity_integralPA15.temporal = sum(singleMeridian_coneDensity(3,:),'omitnan');
singleMeridian_coneDensity_integralPA15.inferior = sum(singleMeridian_coneDensity(4,:),'omitnan');

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

% ------ Visualize HVA and VMA ------ 

titleStr = sprintf('HVA VMA cone density %s - ISETBIO left eye', dataSet);
fH1 = plotHVAandVMA(coneDensity(:,meridians)', eccDeg, titleStr, figureDir, saveFigures);

figure(fH1); set(gca, 'XLim', [0 7]); hold on;
plot(mean([1,6]),hva(coneDataMeridians.singleMeridian), 'ko', 'MarkerSize',10, 'MarkerFaceColor', 'k');
plot(mean([1,6]),vma(coneDataMeridians.singleMeridian), 'ko', 'MarkerSize',10, 'lineWidth',2);

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
    mRGCfDensityInterp(:,aa) = interp1([cardinalMeridianAngles, 360], Watson_mRGCRFDensityPerDeg2(:,aa), angDeg, 'linear');
    totalRGCfDensityInterp(:,aa) = interp1([cardinalMeridianAngles, 360], Watson_totalRGCRFDensityPerDeg2(:,aa), angDeg, 'linear');
end

mRGCfDensity = mRGCfDensityInterp';
totalRGCfDensity = totalRGCfDensityInterp';

% Select mRGC or total RGC in wedge for every meridian
for ii = 1:numel(meridians)
    
    thisWedgeX = wedgeMaskX(:,:,ii);
    thisWedgeY = wedgeMaskY(:,:,ii);
    
    thisSingleMeridianX = singleMeridianX(:,:,ii);
    thisSingleMeridianY = singleMeridianY(:,:,ii);

    wedge_mRGCfWatson(ii,:) = mRGCfDensity(~isnan(thisWedgeX));
    singleMeridian_mRGCfWatson(ii,:) = mRGCfDensity(~isnan(thisSingleMeridianX));
    
    wedge_totalRGCWatson(ii,:) = totalRGCfDensity(~isnan(thisWedgeX));
    singleMeridian_totalRGCWatson(ii,:) = totalRGCfDensity(~isnan(thisSingleMeridianX));
    
end

% Calculate mRGC density integral for +/- 15 deg wedge
wedge_mRGCfWatson_integralPA15.nasal    = sum(wedge_mRGCfWatson(1,:));
wedge_mRGCfWatson_integralPA15.superior = sum(wedge_mRGCfWatson(2,:));
wedge_mRGCfWatson_integralPA15.temporal = sum(wedge_mRGCfWatson(3,:));
wedge_mRGCfWatson_integralPA15.inferior = sum(wedge_mRGCfWatson(4,:));

% Calculate total RGC density integral for +/- 15 deg wedge
wedge_totalRGCWatson_integralPA15.nasal    = sum(wedge_totalRGCWatson(1,:));
wedge_totalRGCWatson_integralPA15.superior = sum(wedge_totalRGCWatson(2,:));
wedge_totalRGCWatson_integralPA15.temporal = sum(wedge_totalRGCWatson(3,:));
wedge_totalRGCWatson_integralPA15.inferior = sum(wedge_totalRGCWatson(4,:));

% Or take the integral of cone density along just the meridia for mRGC
singleMeridian_mRGCfCurcio_integralPA15.nasal    = sum(singleMeridian_mRGCfWatson(1,:));
singleMeridian_mRGCfCurcio_integralPA15.superior = sum(singleMeridian_mRGCfWatson(2,:));
singleMeridian_mRGCfCurcio_integralPA15.temporal = sum(singleMeridian_mRGCfWatson(3,:));
singleMeridian_mRGCfCurcio_integralPA15.inferior = sum(singleMeridian_mRGCfWatson(4,:));

% and total RGC
singleMeridian_totalRGCWatson_integralPA15.nasal    = sum(singleMeridian_totalRGCWatson(1,:));
singleMeridian_totalRGCWatson_integralPA15.superior = sum(singleMeridian_totalRGCWatson(2,:));
singleMeridian_totalRGCWatson_integralPA15.temporal = sum(singleMeridian_totalRGCWatson(3,:));
singleMeridian_totalRGCWatson_integralPA15.inferior = sum(singleMeridian_totalRGCWatson(4,:));


% Concatenate data
mrgcDataMeridians.wedge = [wedge_mRGCfWatson_integralPA15.nasal, ...
                           wedge_mRGCfWatson_integralPA15.superior, ...
                           wedge_mRGCfWatson_integralPA15.temporal, ...
                           wedge_mRGCfWatson_integralPA15.inferior];

mrgcDataMeridians.singleMeridian = [singleMeridian_mRGCfCurcio_integralPA15.nasal, ...
                          singleMeridian_mRGCfCurcio_integralPA15.superior, ...
                          singleMeridian_mRGCfCurcio_integralPA15.temporal, ...
                          singleMeridian_mRGCfCurcio_integralPA15.inferior];
                       
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

% ------ Visualize density and HVA vs VMA ------ 

titleStr = 'HVA VMA mRGCf density Watson 2014 - ISETBIO';
fH2 = plotHVAandVMA(Watson_mRGCRFDensityPerDeg2(1:4,:), eccDeg, titleStr, figureDir, saveFigures);

figure(fH1); h = get(fH2, 'Children');
for ii = [3,2]
    plot(h(2).Children(ii).XData,h(2).Children(ii).YData, ...
        'Color', [0.7 0.7 0.7], ...
        'LineWidth', h(2).Children(ii).LineWidth, ...
        'LineStyle', h(2).Children(ii).LineStyle); hold on
end

figure(fH1); set(gca, 'XLim', [0 7]); hold on;
plot(mean([1,6]),hva(mrgcDataMeridians.singleMeridian),  'o','Color', [0.8 0.8 0.8], 'MarkerSize',10, 'MarkerFaceColor', [0.8 0.8 0.8])
plot(mean([1,6]),vma(mrgcDataMeridians.singleMeridian), 'o', 'Color', [0.8 0.8 0.8], 'MarkerSize',10, 'LineWidth',2)

obj = findobj(gca,'Type','Line');
legend(obj([9,8,4,3]), {'HVA Cones Song et al (2011) - ISETBIO', 'VMA Cones Song et al (2011) - ISETBIO', ...
        'HVA mRGC Watson (2014)', 'VMA mRGC Watson (2014)'}, 'Location', 'SouthEast');


%% -----------------------------------------------------------------
%  --------- mRGC from Curcio et al 1990 rgcDisplacementMap --------
%  -----------------------------------------------------------------

load(fullfile(pfRV1rootPath, 'external', 'data', 'mRFDensityByMeridian.mat'), 'mRFDensityByMeridian', 'regularSupportPosDegVisual','sampleResPolAng')

dt = diff(regularSupportPosDegVisual([1,2]));
eccenLimitsRGC = ((eccenBoundary(1)/dt)+1):1:((eccenBoundary(2)/dt)+1);

mRGC_rgcDisplMap = mRFDensityByMeridian(:, eccenLimitsRGC(1:10:end));

% Interpolate the axes to get 2D density
mRGC_rgcDisplMap(73,:) = mRGC_rgcDisplMap(1,:);

for aa = 1:numel(eccDeg)
    mRGCfDensity_rgcDisplMapInterp(:,aa) = interp1([sampleResPolAng, 360], mRGC_rgcDisplMap(:,aa), angDeg, 'linear');
end

mRGCfDensity_rgcDisplMap = mRGCfDensity_rgcDisplMapInterp';

% Select mRGC or total RGC in wedge for every meridian
for ii = 1:numel(meridians)
    
    thisWedgeX = wedgeMaskX(:,:,ii);
    thisWedgeY = wedgeMaskY(:,:,ii);
    
    thisSingleMeridianX = singleMeridianX(:,:,ii);
    thisSingleMeridianY = singleMeridianY(:,:,ii);

    wedge_mRGCfCurcio(ii,:) = mRGCfDensity_rgcDisplMap(~isnan(thisWedgeX));
    singleMeridian_mRGCfCurcio(ii,:) = mRGCfDensity_rgcDisplMap(~isnan(thisSingleMeridianX));
end

% Calculate mRGC density integral for +/- 15 deg wedge
wedge_mRGCfCurcio_integralPA15.nasal    = sum(wedge_mRGCfCurcio(1,:));
wedge_mRGCfCurcio_integralPA15.superior = sum(wedge_mRGCfCurcio(2,:));
wedge_mRGCfCurcio_integralPA15.temporal = sum(wedge_mRGCfCurcio(3,:));
wedge_mRGCfCurcio_integralPA15.inferior = sum(wedge_mRGCfCurcio(4,:));

% Or take the integral of cone density along just the meridia for mRGC
singleMeridian_mRGCfCurcio_integralPA15.nasal    = sum(singleMeridian_mRGCfCurcio(1,:));
singleMeridian_mRGCfCurcio_integralPA15.superior = sum(singleMeridian_mRGCfCurcio(2,:));
singleMeridian_mRGCfCurcio_integralPA15.temporal = sum(singleMeridian_mRGCfCurcio(3,:));
singleMeridian_mRGCfCurcio_integralPA15.inferior = sum(singleMeridian_mRGCfCurcio(4,:));

% Concatenate data
mrgcDataMeridiansCurcio.wedge = [wedge_mRGCfCurcio_integralPA15.nasal, ...
                           wedge_mRGCfCurcio_integralPA15.superior, ...
                           wedge_mRGCfCurcio_integralPA15.temporal, ...
                           wedge_mRGCfCurcio_integralPA15.inferior];

mrgcDataMeridiansCurcio.singleMeridian = [singleMeridian_mRGCfCurcio_integralPA15.nasal, ...
                          singleMeridian_mRGCfCurcio_integralPA15.superior, ...
                          singleMeridian_mRGCfCurcio_integralPA15.temporal, ...
                          singleMeridian_mRGCfCurcio_integralPA15.inferior];
                       
fprintf('\nmRGCf density Curcio et al (1990), 1-6 deg eccen, +/- 15 deg polar angle wedge around meridia\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(mrgcDataMeridiansCurcio.wedge))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(mrgcDataMeridiansCurcio.wedge))

fprintf('\nmRGCf density Curcio et al (1990), 1-6 deg eccen, only on meridian (no WEDGE)\n')
fprintf('Horizontal-Vertical Asymmetry (Horz / Vert):\t %1.2f%%\n', hva(mrgcDataMeridiansCurcio.singleMeridian))
fprintf('Vertical-Meridian Asymmetry (North / South):  \t %1.2f%%\n', vma(mrgcDataMeridiansCurcio.singleMeridian))

% ------ Visualize density and HVA vs VMA ------ 

titleStr = 'HVA VMA mRGCf density Curcio 1990 - rgcDisplacementMap';
fH3 = plotHVAandVMA(mRGC_rgcDisplMap(ismember(sampleResPolAng,cardinalMeridianAngles),:), eccDeg, titleStr, figureDir, saveFigures);

figure(fH1); h = get(fH3, 'Children');
for ii = [3,2]
    plot(h(2).Children(ii).XData,h(2).Children(ii).YData, ...
        'Color', [0.5 0.5 0.5], ...
        'LineWidth', h(2).Children(ii).LineWidth, ...
        'LineStyle', h(2).Children(ii).LineStyle); hold on
end

figure(fH1); set(gca, 'XLim', [0 7]); hold on;
plot(mean([1,6]),hva(mrgcDataMeridiansCurcio.singleMeridian), 'o','Color', [0.4 0.4 0.4], 'MarkerSize',10, 'MarkerFaceColor', [0.4 0.4 0.4])
plot(mean([1,6]),vma(mrgcDataMeridiansCurcio.singleMeridian), 'o', 'Color', [0.4 0.4 0.4], 'MarkerSize',10, 'LineWidth',2)

obj = findobj(gca,'Type','Line');
legend(obj([13,12, 8,7, 4,3]), {'HVA Cones Song et al (2011) - ISETBIO', 'VMA Cones Song et al (2011) - ISETBIO', ...
        'HVA mRGC Watson (2014) - ISETBIO', 'VMA mRGC Watson (2014) - ISETBIO', ...
        'HVA mRGC Curcio (1990) - rgcDisplMap', 'VMA mRGC Curcio (1990) - rgcDisplMap'}, 'Location', 'SouthEast');
title('Cone and mRGC density asymmetries')
savefig(fullfile(figureDir, 'Cone_RGC_HVA_VMA_eccen1-6_VSS2020_Benson'))
print(fullfile(figureDir, 'Cone_RGC_HVA_VMA_eccen1-6_VSS2020_Benson'), '-dpdf', '-fillpage')


% ---- Bar graph of VMA and VMA -----

ttl = 'Density integral of 15^o wedges, 1-6^o eccen';
hvas = [hva(coneDataMeridians.wedge), hva(mrgcDataMeridiansCurcio.wedge)];
vmas = [vma(coneDataMeridians.wedge), vma(mrgcDataMeridiansCurcio.wedge)];

pos = [1, 1.1, 1.3, 1.4];

figure(100); clf; set(gcf, 'Position', [230   273   370   437], 'Color', 'w');
bar(pos(1), hvas(1),  0.1, 'EdgeColor','k','facecolor','k', 'FaceAlpha', 0.8, 'EdgeAlpha', 0.8, 'LineWidth',3); hold on
bar(pos(2), vmas(1),  0.1, 'EdgeColor','k','facecolor','none', 'FaceAlpha', 0.8, 'EdgeAlpha', 0.8,'LineWidth',3); hold on
bar(pos(3), hvas(2),  0.1, 'EdgeColor',[0.7 0.7 0.7],'facecolor',[0.7 0.7 0.7], 'FaceAlpha', 0.8, 'EdgeAlpha', 0.8,'LineWidth',3); hold on
bar(pos(4), vmas(2),  0.1, 'EdgeColor',[0.7 0.7 0.7],'facecolor','none', 'FaceAlpha', 0.8,'EdgeAlpha', 0.8,'LineWidth',3); hold on
title(ttl)
legend({'HVA', 'VMA'}, 'Location', 'SouthEast'); legend boxoff;
set(gca,'xlim',[0.9,1.5],'ylim',[-10,27], 'FontSize', 15);
set(gca,'XTick', [1.05,1.35], 'XTickLabel',{'Cones', 'mRGC'}, 'TickDir', 'out');
ylabel('Asymmetry %'); box off;

savefig(fullfile(figureDir, 'Cone_RGC_HVA_VMA_eccen1-6_VSS2020_Benson_summary'))
print(fullfile(figureDir, 'Cone_RGC_HVA_VMA_eccen1-6_VSS2020_Benson_summary'), '-dpdf', '-fillpage')

