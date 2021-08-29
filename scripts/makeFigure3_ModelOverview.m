% makeFigure3_ModelOverview


% Set paths
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensitynophaseshiftlonly500';
figPth     = fullfile(baseFolder, 'figures', 'Figure3_modeloverview');
saveFigs   = true; % Save figures?

% dataPth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/';
expstr     = 'contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat';

% load the simulated photocurrent data
load(fullfile(figPth,sprintf('Mn_OGconeOutputs_%s', expstr)));

% Get size
m2deg = 1./expParams.deg2m;
xPosInDeg = cMosaic.size * m2deg;

% Get colors
cmap = [1 0 0; 0 1 0; 0 0 1];


%% 1. RADIANCE

% load stimulus data
load(fullfile(baseFolder,'figures','Figure3_modeloverview','stimulus_100.mat'));

figure(1); clf; set(gcf, 'Color','w', 'Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Radiance');

radiance = scenes{2}.data.photons;
midpoint = ceil(size(radiance,1)/2);

% Plot scene radiance
imagesc(sum(radiance,3));
% hold on; plot([0 size(radiance,1)], [midpoint midpoint], 'k:');
set(gca,'CLim', [0, 2.3].*10^17, 'TickDir','out', 'FontSize',15);
colormap gray; colorbar; axis image; box off;
title('Scene radiance summed over nm','Fontsize',12)
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    print(fullfile(figPth, 'Fig3_radiance_2d'),'-depsc')
end


%% 2. IRRADIANCE
figure(2); clf; set(gcf, 'Color','w', 'Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Irradiance');

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
colormap(cInt); axis image; cb = colorbar; box off;
set(gca,'CLim', [0 1], 'TickDir','out', 'FontSize',15); cb.Ticks = linspace(0, 1, 4) ; %Create 8 ticks from zero to 1
cb.TickLabels = num2cell([ylI(1), ylI(2)*(1/3),  ylI(2)*(2/3), ylI(2)]);
title('Retinal irradiance','Fontsize',12);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')


if saveFigs
    print(fullfile(figPth, 'Fig3_modeloverview_irradiance_2d'),'-dpng')
end


%% Plot mosaic pattern
figure(3); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Mosaic');
title('Cone Mosaic')
[X,Y] = meshgrid(1:cMosaic.rows,1:cMosaic.cols);
plot(X,Y,'ko')
set(gca,'CLim', [2 4]); axis square; box off;
colormap(cmap)
midpoint = ceil(cMosaic.rows/2);
quarterpoint = ceil(cMosaic.rows/4);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')


if saveFigs
    hgexport(3,fullfile(figPth, 'Fig3_modeloverview_coneMosaic.eps'))
end

%% Load absorptions, current and RGC responses for 4.5 deg eccen 100% contrast - L cone only
load(fullfile(figPth,'rgcResponses_latenoiselevel1.0_withDownsampling_c32_eccen4.50.mat'))

%% Absorptions
midpoint = ceil(size(x.absorptions,1)/2);
figure(4); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Absorptions');
imagesc(squeeze(x.absorptions(:,:,1)));
axis square; title('Cone absorptions');
colormap gray; colorbar; box off;
set(gca,'CLim',[0 500], 'TickDir','out');
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(4,fullfile(figPth, 'Fig3_modeloverview_coneAbsorptions.eps'))
end

%% Current
figure(5); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Current');
imagesc(squeeze(x.current(:,:,1)));
colormap gray; colorbar; box off;
axis square; title('Current');
set(gca,'CLim',[-35 0], 'TickDir','out');
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(5,fullfile(figPth, 'Fig3_modeloverview_coneCurrent.eps'))
end

%% Current filtered
figure(6);  clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Current Filtered by RGC DoG');
imagesc(squeeze(x.Filtered(:,:,1)));
colormap gray; colorbar; box off;
axis square; title('Filtered');
set(gca,'CLim',[-35 0]);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(6,fullfile(figPth, 'Fig3_modeloverview_currentRGCFiltered.eps'))
end

%% Current filtered + late Noise
figure(7);  clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - RGC Response (w/ Late Noise)');
imagesc(squeeze(x.LateNoise(:,:,1)));
colormap gray; colorbar; box off;
axis square; title('Late Noise');
set(gca,'CLim',[-35 0]);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(7,fullfile(figPth, 'Fig3_modeloverview_RGCFilteredLateNoise.eps'))
end

%% RGC response mRGC:cone = 2:1
figure(8); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - RGC Downsampled 2:1');
imagesc(squeeze(x.DownSampled1(:,:,1)));
colormap gray; colorbar;
axis square; title('DownSampled1');
set(gca,'CLim',[-35 0]);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(8,fullfile(figPth, 'Fig3_modeloverview_RGCFilteredLateNoiseDownsampled1.eps'))
end

%% RGC response mRGC:cone = 1:0.08
figure(9); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - RGC Downsampled 1:0.08');
imagesc(squeeze(x.DownSampled5(:,:,1)));
colormap gray; colorbar;
axis square; title('DownSampled5');
set(gca,'CLim',[-35 0]);
xlabel('x position (deg)')
ylabel('y position (deg)')
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(9,fullfile(figPth, 'Fig3_modeloverview_RGCFilteredLateNoiseDownsampled5.eps'))
end


%% DoG filter

cRows = size(x.absorptions,1);
cCols = size(x.absorptions,2);
colors = parula(6);

figure(10); hold all; set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - FTT RGC DoG filter'); clf; hold all;
box off; axis square;
xlabel('Spatial Frequency (cycles/deg)');
ylabel('Normalized amplitude (a.u.)');
set(gca,'TickDir','out','FontSize',15)

for cone2RGCRatio = 1:5

    sigma.center   = cone2RGCRatio*rgcParams.DoG.kc; % Center Gauss RGC
    sigma.surround = cone2RGCRatio*rgcParams.DoG.ks; % Surround Gauss RGC
    sigma.ratio    = sigma.surround/sigma.center;              % ratio center surround
    vol.ratio      = rgcParams.DoG.ws/rgcParams.DoG.wc;        % ratio of the surround volume to the center volume

    % Create RGC grid by resampling with cone2RGC ratio.
    rowIndices = 1:cone2RGCRatio:cRows;
    colIndices = 1:cone2RGCRatio:cCols;
    rgcarray = zeros(cRows, cCols);
    rgcarray(rowIndices, colIndices) = 1;

    % Create DoG filter
    [DoGfilter,xx,yy] = makedog2d(31,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);
    DoGfilter = DoGfilter / sum(DoGfilter(:));
    midpoint = ceil(size(DoGfilter,1)/2);
    fftDoG = abs(fft(DoGfilter(midpoint,:)));
    xax = linspace(-midpoint,midpoint,length(fftDoG));
    % Plot if
    figure(10); 
    plot(xax, fftshift(fftDoG), 'color', colors(cone2RGCRatio,:), 'lineWidth',3); hold on;
    plot([4,4], [-0.1 3], 'r:', 'lineWidth',2);
%     plot([0,length(fftDoG)], [0 0], 'k');
    xlim([0 midpoint]); ylim([0, 2.5])
    
    
    if ismember(cone2RGCRatio,[1,5])
        % Plot 1D Filter
        figure(10+cone2RGCRatio); clf; set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC DoG filter'); clf; hold all;
        sz = size(DoGfilter);
        midpoint = ceil(sz(1)/2);
        plot(DoGfilter(midpoint,:), 'color', 'k', 'LineWidth',2); hold all;
        plot([0,sz(1)], [0 0], 'k'); hold on;

        set(gca,'YLim', [-0.05,2], 'XLim', [0, sz(1)], 'TickDir', 'out', 'FontSize', 15, 'XTick', [0, midpoint, sz(1)], ...
            'XTickLabel', {'-1', '0', '1'});
        xlabel('x position (deg)')
        ylabel('Modulation (a.u.)')
        title('DoG filter')
        axis square; box off;

        if saveFigs
            hgexport(10+cone2RGCRatio, fullfile(figPth, sprintf('Fig3_modeloverview_DoGfilter1D_ratio%d',cone2RGCRatio)))
        end
        
        if cone2RGCRatio==5
            set(gca,'YLim', [-0.02,0.15])

            if saveFigs
                hgexport(10+cone2RGCRatio, fullfile(figPth, sprintf('Fig3_modeloverview_DoGfilter1D_size%d_yscaled',cone2RGCRatio)))
            end
        end
        
        % Plot 2D Filter small
        figure(11+cone2RGCRatio); clf; set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC DoG filter'); clf; hold all;

        th = 0:pi/100:(2*pi);
        Xcenter_small = sigma.center * cos(th);
        Ycenter_small = sigma.center * sin(th);
        Xsurround_small = sigma.surround * cos(th);
        Ysurround_small = sigma.surround * sin(th);
        plot(Xcenter_small,Ycenter_small, 'color', 'r', 'LineWidth',2); hold all;
        plot(Xsurround_small,Ysurround_small, 'color', 'g', 'LineWidth',2); hold all;
        plot([-midpoint,midpoint], [0, 0], 'k',  [0 0], [-midpoint,midpoint], 'k'); hold on;
        xlim([-midpoint midpoint]); ylim([-midpoint midpoint])

        set(gca, 'TickDir', 'out', 'FontSize', 15, ...
            'XTick', [-midpoint, 0, midpoint], ...
            'XTickLabel', {'-1', '0', '1'}, ...
            'YTick', [-midpoint, 0, midpoint], ...
            'YTickLabel', {'-1', '0', '1'});
        xlabel('x position (deg)')
        xlabel('y position (deg)')

        title('DoG filter')
        axis square; box off;

        if saveFigs
            hgexport(11+cone2RGCRatio, fullfile(figPth, sprintf('Fig3_modeloverview_DoGfilter2D_size%d',cone2RGCRatio)))
        end



    end
end

if saveFigs
    hgexport(10, fullfile(figPth, sprintf('Fig3_modeloverview_FFTDoGfilter1D')))
end