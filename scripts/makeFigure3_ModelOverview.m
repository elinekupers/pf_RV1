% makeFigure3_ModelOverview

% Set paths
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
expName    = 'conedensitynophaseshiftlonly500';
figPth     = fullfile(baseFolder, 'figures', 'Figure3_modeloverview');
saveFigs   = false; % Save figures?

% dataPth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/';
expstr     = 'contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-1.00.00.0.mat';

% load the simulated photocurrent data
load(fullfile(figPth,sprintf('OGconeOutputs_%s', expstr)));
load(fullfile(figPth,sprintf('current_OGconeOutputs_%s', expstr)));

% Load absorptions, current and RGC responses for 4.5 deg eccen 100% contrast - L cone only
rgc_eccen4_5_c100 = load(fullfile(figPth,'rgcResponses_latenoiselevel1.0_withDownsampling_c32_eccen4.50.mat'));
rgc_eccen10_c100 = load(fullfile(figPth,'rgcResponses_latenoiselevel1.0_withDownsampling_c32_eccen10.00.mat'));
rgc_eccen4_5_c10 = load(fullfile(figPth,'rgcResponses_latenoiselevel1.0_withDownsampling_c22_eccen4.50.mat'));
rgc_eccen40_c100 = load(fullfile(figPth,'rgcResponses_latenoiselevel1.0_withDownsampling_c32_eccen40.00.mat'));

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
midpnt = ceil(size(radiance,1)/2);

% Plot scene radiance
imagesc(sum(radiance,3));
% hold on; plot([0 size(radiance,1)], [midpoint midpoint], 'k:');
set(gca,'CLim', [0, 2.3].*10^17, 'TickDir','out', 'FontSize',15);
colormap gray; colorbar; axis image; box off;
title('Scene radiance summed over nm','Fontsize',12)
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
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
midpnt   = ceil(size(irradiance,1)/2);
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
    'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
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
midpnt = ceil(cMosaic.rows/2);
quarterpoint = ceil(cMosaic.rows/4);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')


if saveFigs
    hgexport(3,fullfile(figPth, 'Fig3_modeloverview_coneMosaic.eps'))
end

%% Single time points for Absorptions
midpnt = ceil(size(absorptions,2)/2);

for t = [1:3]
    figure(4); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Absorptions');
    imagesc(squeeze(absorptions(1,:,:,t,1)));
    axis square; title(sprintf('Cone absorptions t=%d',t));
    colormap gray; colorbar; box off;
    set(gca,'CLim',[0 500], 'TickDir','out');
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
        'XTickLabel',{'-0.5', '0', '0.5'}, ...
        'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
        'YTickLabel',{'-0.5', '0', '0.5'});
    xlabel('x position (deg)')
    ylabel('y position (deg)')

    if saveFigs
        hgexport(4,fullfile(figPth, sprintf( 'Fig3_modeloverview_coneAbsorptions_t%d.eps',t)))
    end

end

%% Single time points for Current
midpnt = ceil(size(current,2)/2);

for t = [1:3]
    figure(104); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Current');
    imagesc(squeeze(current(1,:,:,t+50,1)));
    axis square; title(sprintf('Cone current t=%d',t));
    colormap gray; colorbar; box off;
    set(gca,'CLim',[-40 0], 'TickDir','out');
    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
        'XTickLabel',{'-0.5', '0', '0.5'}, ...
        'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
        'YTickLabel',{'-0.5', '0', '0.5'});
    xlabel('x position (deg)')
    ylabel('y position (deg)')

    if saveFigs
        hgexport(104,fullfile(figPth, sprintf( 'Fig3_modeloverview_coneCurrent_t%d.eps',t)))
    end

end

%% Current filters

figure(105); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - Current');
plot(cMosaic.interpFilterTimeAxis,interpFilters(:,1), 'r', 'lineWidth',2); hold on;
plot([0 max(cMosaic.interpFilterTimeAxis)],[0 0], 'k')
axis square; title(sprintf('L-cone filter for current'));
box off;
set(gca, 'TickDir', 'out', 'FontSize', 15)
xlabel('time (s)')
ylabel('Modulation (a.u.)')

if saveFigs
    hgexport(105,fullfile(figPth, sprintf( 'Fig3_modeloverview_currentFilters.eps')))
end

%% RGC response mRGC:cone = 2:1
figure(8); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - RGC Downsampled 2:1');
imagesc(squeeze(rgc_eccen4_5_c100.x.DownSampled1(:,:,1)));
colormap gray; colorbar;
axis square; title('DownSampled1');
set(gca,'CLim',[-35 0]);
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(8,fullfile(figPth, 'Fig3_modeloverview_RGCFilteredLateNoiseDownsampled1.eps'))
end

%% RGC response mRGC:cone = 1:0.08
figure(9); clf; set(gcf, 'Color', 'w','Position', [459,41,875,737], 'NumberTitle','off', 'Name', 'Fig 3 - RGC Downsampled 1:0.08');
imagesc(squeeze(rgc_eccen4_5_c100.x.DownSampled5(:,:,1)));
colormap gray; colorbar;
axis square; title('DownSampled5');
set(gca,'CLim',[-35 0]);
xlabel('x position (deg)')
ylabel('y position (deg)')
set(gca, 'TickDir', 'out', 'FontSize', 15, ...
    'XTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'XTickLabel',{'-0.5', '0', '0.5'}, ...
    'YTick', [midpnt*0.5, midpnt, midpnt*1.5], ...
    'YTickLabel',{'-0.5', '0', '0.5'});
xlabel('x position (deg)')
ylabel('y position (deg)')

if saveFigs
    hgexport(9,fullfile(figPth, 'Fig3_modeloverview_RGCFilteredLateNoiseDownsampled5.eps'))
end


%% DoG filter 4.5 deg eccen
rgcParams      = rgc_eccen4_5_c100.rgcParams;
sigma.center   = rgcParams.DoG.kc; % Center Gauss RGC
sigma.surround = rgcParams.DoG.ks; % Surround Gauss RGC
sigma.ratio    = sigma.surround/sigma.center;              % ratio center surround
vol.ratio      = rgcParams.DoG.ws/rgcParams.DoG.wc;        % ratio of the surround volume to the center volume

% Create DoG filter
supportDoG = 31;
[DoGfilter,xx,yy] = makedog2d(supportDoG,[],[],sigma.center,sigma.ratio,vol.ratio,[],[]);
DoGfilter = DoGfilter / sum(DoGfilter(:));
midpnt = ceil(size(DoGfilter,1)/2);
fftDoG = fftshift(abs(fft(DoGfilter(midpnt,:))));

conesPerRow = [size(rgc_eccen4_5_c100.x.absorptions,1), ...
               size(rgc_eccen10_c100.x.absorptions,1), ...
               size(rgc_eccen40_c100.x.absorptions,1)];

for cRows = conesPerRow

    % Compute cones / deg
    x              = -15:15;                % DoG support  % 31 cones per image
    fs             = 0:length(x)-1;         % cycles per image
    conesperdegree = cRows / cparams.cmFOV;  % 39.5 cones per deg
    degperimage    = 1/conesperdegree * supportDoG; % 1/39 (deg per cone) * 31 (cones / per image);
    fs             = fs / degperimage;      % cycles / deg

    % Plot 1D Filter
    figure(11); clf; 
    set(gcf, 'Position', [159,382,1398,406], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - RGC DoG filter'); clf; hold all;
    subplot(131)
    plot(DoGfilter(midpnt,:), 'color', 'k', 'LineWidth',2); hold all;
    plot([1,size(DoGfilter,1)], [0 0], 'k'); hold on;
    set(gca,'YLim', [-0.1,2], 'XLim', (size(DoGfilter,1)/2)+[-conesperdegree/4,conesperdegree/4], ...
        'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [(size(DoGfilter,1)/2)-conesperdegree/4,ceil((size(DoGfilter,1)/2)),(size(DoGfilter,1)/2)+conesperdegree/4], ...
        'XTickLabel', [-1, 0, 1]);
    xlabel('x position (deg)')
    ylabel('Modulation (a.u.)')
    title('DoG 1D visual space')
    axis square; box off;

    %% Plot 2D Filter small

    th = 0:pi/100:(2*pi);
    Xcenter   = 2*sigma.center * cos(th);
    Ycenter   = 2*sigma.center * sin(th);
    Xsurround = 2*sigma.surround * cos(th);
    Ysurround = 2*sigma.surround * sin(th);
    subplot(132); cla;
    for ii = ceil(-conesperdegree/2):1:ceil(conesperdegree/2)
        plot([ii, ii],[-conesperdegree/2, conesperdegree/2],'color',[0.7 0.7 0.7]); hold all;
        plot([-conesperdegree/2, conesperdegree/2],[ii, ii],'color',[0.7 0.7 0.7]); hold all;
    end
    plot([-conesperdegree/2 conesperdegree/2], [0, 0], 'k', ...
        [0,0],[-conesperdegree/2 conesperdegree/2], 'k');
    plot(Xcenter,Ycenter, 'color', 'r', 'LineWidth',2); hold all;
    plot(Xsurround,Ysurround, 'color', 'b', 'LineWidth',2); 
    xlim([-conesperdegree/2 conesperdegree/2]); ylim([-conesperdegree/2 conesperdegree/2])

    set(gca, 'TickDir', 'out', 'FontSize', 15, ...
        'XTick', [-conesperdegree/2, 0, conesperdegree/2], ...
        'XTickLabel', [-0.5, 0, 0.5], ...
        'YTick', [-conesperdegree/2, 0, conesperdegree/2], ...
        'YTickLabel', [-0.5, 0, 0.5]);
    xlabel('x position (deg)')
    xlabel('y position (deg)')

    title('DoG 2D visual space')
    axis square; box off;


    %% Plot 1D Filter
    subplot(133); 
    plot(fs(1:15),fftDoG(1:15), 'color', 'k', 'LineWidth',2); hold all;
    plot([4,4],[0,2],'r:','lineWidth',2)
    set(gca,'YLim', [0, 2], 'XLim', [0 18], ...
        'TickDir', 'out', 'FontSize', 15);
    xlabel('Spatial Frequency (cycles/deg)')
    ylabel('Amplitude (a.u.)')
    title('DoG 1D FFT')
    axis square; box off;


    if saveFigs
        hgexport(11, fullfile(figPth,sprintf('Fig3_modeloverview_DoGfilter1D_2D_fft_coneRows%d',cRows)))
    end
end


%% SVM weights

% label = [ones(1000, 1); -ones(1000, 1)];
% data = permute(rgc_eccen4_5_c100.x.DownSampled1, [4,1,2,3]);
% data = reshape(data, 2000, []);
% 
% % Fit the SVM model.
% cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
% 
% d = 0.02;
% [x1Grid,x2Grid] = meshgrid(min(data(1,:)):d:max(data(1,:)), ...
%                     min(data(1001,:)):d:max(data(1001,:)));
% xGrid = [x1Grid(:),x2Grid(:)];
% 
% figure(13); clf; set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - SVM betas');    
% h = nan(3,1); % Preallocation
% cw = reshape(data(1,:),[79 79]);
% ccw = reshape(data(1001,:),[79 79]);
% grp = [ones(size(cw(:,39)));-ones(size(cw(:,39)))];
% hold all;
% h(1:2) = gscatter([cw(39,:),ccw(39,:)],[cw(:,39);ccw(:,39)],grp,'rg','+*');
% h(3) = plot(cdata(cvmdl.IsSupportVector,1),...
%     cdata(cvmdl.IsSupportVector,2),'ko');
% contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
% legend(h,{'-1','+1','Support Vectors'},'Location','Southeast');
% axis equal
% hold off
% 
% 
% figure(13); clf; set(gcf, 'Position', [786,135,560,420], 'color', 'w', 'NumberTitle','off', 'Name', 'Fig 3 - SVM betas');    
% hold all;
% cw = reshape(data(1,:),[79,79]);
% ccw = reshape(data(1001,:),[79,79]);
% gscatter([cw(39,:),ccw(39,:)],[cw(:,39);ccw(:,39)]',[ones(1,79),-ones(1,79)]); 

% plot(cvmdl.CrossValidatedModel.
box off; set(gca, 'TickDir', 'out', 'FontSize',12); 
colormap gray; axis image; colorbar;
  set(gca, 'TickDir', 'out', 'FontSize', 15, ...
            'XTick', [1, 40, 79], ...
            'XTickLabel', {'-1', '0', '1'}, ...
            'YTick', [1, 40, 79], ...
            'YTickLabel', {'-1', '0', '1'});
set(gca,'CLim', 10e-4*[-0.3624,0.3624]);
    print(fullfile(figPth, 'classifierWeights_averagedAcrossTime.eps'),'-deps')
    
    