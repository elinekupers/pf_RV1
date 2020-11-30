% makeFigure3_ModelOverview

baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName = 'default';
ratios = 1:5;
labelsRatio = sprintfc('cone2RGC = %1.1f:1', 2./ratios);
colors = parula(6);

% Save figures?
saveFigs = false;

saveFigPth = fullfile(baseFolder, 'figures', 'Figure3', 'Figure3_modeloverview');
if ~exist(saveFigPth); mkdir(saveFigPth); end

%% Plot the cone mosaic separate because of the different colormap

cmap = [1 0 0; 0 1 0; 0 0 1];
varySatValues = 1-linspace(0.1,1,5);
trialColors = varysat(cmap,varySatValues);
contrastColors = [0 0 0; 241 101 33]./255;

figPath = fullfile(ogRootPath, 'figs', 'overviewFigure5');
load(fullfile(figPath,'cMosaicPlotted'))

% Get size
m2deg = 10/3;
xPosInDeg = cMosaic.size .* m2deg;

figure(100); clf; set(gcf, 'Color', 'w','Position', [239, 369, 460, 409], 'NumberTitle','off', 'Name', 'Fig 3 - Cone Mosaic');

% Plot mosaic pattern
title('Cone Mosaic')
imagesc(cMosaic.pattern)
set(gca,'CLim', [2 4]); axis image; axis off
colormap(cmap)

if saveFigs
    hgexport(100,fullfile(saveFigPth, 'Fig3_modeloverview_coneMosaic.eps'))
end


%% Set up other figures

% Select a trial for cone absorptions
trialNum = 4;
stimulusNum = 4;

timeWindowA = [1 2 3];  % select time points for absorptions (ms)
climsA = [0 220];  % color bar / ylims for absorptions (photon count)
ylA = [0 500];  % y limit for time series (absorptions, photon count)

% File names
fNames = {'absorptionsPlotted', 'stimulus'};
contrast = 100;

trialsToPlot = ceil(100*rand(5,1));

% load data
for ii = 1:length(fNames)
    load(fullfile(figPath, [fNames{ii} sprintf('_%d.mat',contrast)]));
end

%% 1. RADIANCE
figure(1); clf; set(1, 'Color','w', 'Position', [239, 369, 460, 409], 'NumberTitle','off', 'Name', 'Fig 3 - Radiance');

radiance = scenes{2}.data.photons;
midpoint = ceil(size(radiance,1)/2);

% Plot scene radiance
imagesc(sum(radiance,3));
hold on; plot([0 size(radiance,1)], [midpoint midpoint], 'k:');
set(gca,'CLim', [0, 2.3].*10^17);
colormap gray; colorbar; axis image; axis off
title('Scene radiance summed over nm','Fontsize',12)

if saveFigs
    print(fullfile(saveFigPth, 'Fig3_radiance_2d'),'-depsc')
end


%% 2. IRRADIANCE
figure(2); clf; set(2, 'Color','w', 'Position', [239, 369, 460, 409], 'NumberTitle','off', 'Name', 'Fig 3 - Irradiance');

% Plot the stimulus after optics
irradiance = OG(2).oiModulated.data.photons;
midpoint   = ceil(size(irradiance,1)/2);
gamma      = 1;
wList      = (400:10:700);

XYZ = ieXYZFromPhotons(irradiance, wList);
XYZ_normalized = XYZ/max(XYZ(:));
rgbIrradianceData = xyz2srgb(XYZ_normalized);

rgbColorMapMin = min(reshape(rgbIrradianceData, [size(rgbIrradianceData,1)*size(rgbIrradianceData,2), 3]));
rgbColorMapMax = max(reshape(rgbIrradianceData, [size(rgbIrradianceData,1)*size(rgbIrradianceData,2), 3]));
cInt = [linspace(rgbColorMapMin(1), rgbColorMapMax(1), 64); ...
    linspace(rgbColorMapMin(2), rgbColorMapMax(2), 64); ...
    linspace(rgbColorMapMin(3), rgbColorMapMax(3), 64)]';

irradianceSumAllWV = sum(irradiance,3);
ylI = [min(irradianceSumAllWV(:)), max(irradianceSumAllWV(:))];

imagesc(rgbIrradianceData);
colormap(cInt); axis image; axis off; cb = colorbar;
set(gca,'CLim', [0 1]); cb.Ticks = linspace(0, 1, 4) ; %Create 8 ticks from zero to 1
cb.TickLabels = num2cell([ylI(1), ylI(2)*(1/3),  ylI(2)*(2/3), ylI(2)]);
title('Retinal irradiance','Fontsize',12)

if saveFigs
    print(fullfile(saveFigPth, 'Fig3_modeloverview_irradiance_2d'),'-dpng')
end


%% 3. ABSORPTIONS

figure(3); clf; set(3, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 3 - Absorptions');

% get average absorptions across time. Integration time was 2 ms, thus
% recalculate to absorption per 1 ms
theseAbsorptions = squeeze(absorptions(trialNum,:,:,timeWindowA,stimulusNum))./2;

midpoint = ceil(size(theseAbsorptions,1)/2);

subplotIdx = [1:3];

for ii = subplotIdx
    subplot(1,3,ii); hold on;
    imagesc(theseAbsorptions(:,:,mod(ii,4)));
    colormap gray; axis image; axis off;
    set(gca,'CLim',climsA); colorbar;
    title(sprintf('Absorptions (photons) at t=%d',mod(ii-1,3)+1),'Fontsize',12)
    
    colormap gray; axis image; axis off;
    set(gca,'CLim',climsA);
end


if saveFigs
    print(fullfile(saveFigPth, 'Fig3_modeloverview_absorptions_2d'),'-dpng')
end


%% All DoG filters in one figure
fH4 = figure(4); set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC DoG filters'); clf; hold all;

for r = ratios
    fNameFilter = fullfile(baseFolder, 'data', expName, 'rgc', sprintf('rgcDoGFilter_Cones2RGC%d_absorptionrate.mat',r));
    load(fullfile(fNameFilter));
    
    sz = size(DoGfilter);
    midpoint = ceil(sz(1)/2);
    plot(DoGfilter(midpoint,:), 'color', colors(r,:), 'LineWidth',2); hold all;
end

% Plot x-axis
plot([0,sz(1)], [0 0], 'k'); hold on;

set(gca,'YLim', [-0.05,1], 'XLim', [0, sz(1)], 'TickDir', 'out', 'FontSize', 15, 'XTick', [0, midpoint, sz(1)], ...
    'XTickLabel', {'-1', '0', '1'});
legend(labelsRatio, 'Location','Best'); legend boxoff;
xlabel('x position (deg)')
ylabel('Normalized modulation (a.u.)')
title('DoG filters')
axis square; box off;

if saveFigs
    hgexport(fH4, fullfile(saveFigPth, 'DoGs_allInOne'))
end

%% 2D RGC RESPONSE FILTERED + SUBSAMPLE
fH5 = figure(5); set(gcf, 'Position', [33,42,1639,763], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC response'); clf; hold all;
r = ratios(1);

fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', sprintf('ratio%d',r), sprintf('rgcResponse_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
load(fullfile(fNameFiltSubResponse));

% Use selected time points only
filteredSubsampled = squeeze(mean(rgcResponse(:,:,:,[1:3],1),1))./2;
midpoint = ceil(size(filteredSubsampled,1)/2);

for ii = 1:3
    subplot(1,3,ii);
    imagesc(squeeze(filteredSubsampled(:,:,ii))); colormap gray; hold on;
    colorbar;
    xlabel('x position (deg)')
    ylabel('y position (deg)')
    title({labelsRatio{r} ['time sample ' ii]});
    set(gca, ...
        'CLim', [0 max(filteredSubsampled(:))], 'TickDir', 'out', 'FontSize', 15)
    axis square; box off;
end


hgexport(fH5, fullfile(saveFigPth, '2D_rgcResponseSubsampled'))



%% 1D and 2D RGC RESPONSE FFT

fH6 = figure(6); set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC 2D FFT'); clf; hold all;

r = ratios(1);

fNameFiltSubResponse = fullfile(baseFolder, 'data', expName, 'rgc', 'figure4', sprintf('ratio%d',r), sprintf('rgcResponse_Cones2RGC%d_contrast1.0000_eccen4.50_absorptionrate.mat',r));
load(fullfile(fNameFiltSubResponse));

% Use selected time points only
filteredSubsampled = squeeze(mean(rgcResponse(:,:,:,1,1),1))./2;
filteredSubsampled_amps_2D  = abs(fft2(filteredSubsampled));
filteredSubsampled_amps_2D(1,1) = NaN;
filteredSubsampled_amps_2D = fftshift(filteredSubsampled_amps_2D);

midpoint = 1+ceil(size(filteredSubsampled_amps_2D,1)/2);
quarterpoint_rsp = midpoint/2;

imagesc(filteredSubsampled_amps_2D); colormap gray;
colorbar;
set(gca, 'CLim', [-1 1].* max(filteredSubsampled_amps_2D(:)))
xlabel('sf (cycles/deg)')
ylabel('sf (cycles/deg)')
title(labelsRatio{r})

set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{num2str(-quarterpoint_rsp/rgcParams.fov), '0', num2str(quarterpoint_rsp/rgcParams.fov)});
axis square; box off;


hgexport(fH6, fullfile(saveFigPth, '2D_rgcResponseFFT'))



